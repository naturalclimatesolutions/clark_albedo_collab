% NCS - Global Albedo - TNC-Clark Collaboration - Import base data
%
% This routine is used to import base data for the Global Albedo calculation, namely MODIS land
%   cover, snow cover data, radiation data, radiative kernels, and albedo atlas data. All data is
%   read, rescaled to 0.05x0.05 lat/lon scale if needed and atlas values are gap-filled as well. It
%   also calculated the land area in each pixel (in WGS84).
% 
% References:
%   Land-cover:     MODIS/Terra – Land Cover - https://lpdaac.usgs.gov/products/mcd12c1v006/
%   Snow-cover:     MODIS/Terra – Snow cover from NSIDC - https://nsidc.org/data/MOD10CM/versions/6
%   Solar radiation:NCAR-NCEP Reanalysis for solar radiation
%   Albedo:         Gao, F., T. He, Z. Wang, B. Ghimire, Y. Shuai, J. Masek, C. Schaaf, 
%                       C. Williams (2014), Multiscale climatological albedo look-up maps derived
%                       from moderate resolution imaging spectroradiometer BRDF/albedo products, 
%                       Journal of Applied Remote Sensing, 8(1), 083532.
%   Kernels:
%                   CACK v1.0 – Bright R. M. and O’Halloran, T. L.: CACKv1.0,
%                       available at: https://doi.org/10.6073/pasta/d77b84b11be99ed4d5376d77fe0043d8
%                   CAM3 – based on [Shell et al., 2008], people.oregonstate.edu/~shellk/kernel.html
%                       (Note 9.29.21 - dataset seems to have been removed)
%                   CAM5 – Pendergrass, A. G.: CAM5 Radiative Kernels
%                       https://www.earthsystemgrid.org/dataset/ucar.cgd.ccsm4.cam5-kernels.html
%                   ECHAM6 - Block, K. and Mauritsen, T.: ECHAM6 CTRL kernel,
%                       https://swiftbrowser.dkrz.de/public/
%                           dkrz_0c07783a-0bdc-4d5e-9f3b-c1b86fac060d/Radiative_kernels/
%                   HadGEM2 – Smith, C. J.: HadGEM2 radiative kernels, edited, University of Leeds.
%                       https://archive.researchdata.leeds.ac.uk/382/
%                   HadGEM3 – Smith, C. J.: HadGEM3-GA7.1 radiative kernels
%                       https://doi.org/10.5281/zenodo.3594673
%   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define storage arrays
% *********************
snowcover_all = zeros(nlat05,nlon05,nmonths);
direct_solar = zeros(nlat1dg,nlon1dg,nmonths);
diffuse_solar = zeros(nlat1dg,nlon1dg,nmonths);
solar_direct_hires  = zeros(nlat,nlon,nmonths);
solar_diffuse_hires = zeros(nlat,nlon,nmonths);
frac_diff_hires     = zeros(nlat,nlon,nmonths);
albkern_hires = zeros(nlat,nlon,nmonths);
MisSnow = false(nlat,nlon,nmonths);     % Used to record where I gap-filled the snow data

% Import land cover data
% **********************
% Note: IGBP classification is given as "land cover type 1"
%       I import the land cover proportion for each land cover-type ("Land_Cover_Type_1_Percent") to
%       calculate albedo-induced radiative forcing on every land pixel that has some land cover.
% Note: MODIS/Terra grid is in 0.05 x 0.05 degree grid, with upperleft corner (90,-180) lat/lon
LCP = hdfread(landcovermap,"Land_Cover_Type_1_Percent");
landmask = true(nlat,nlon);
i = biomes_igbp =="WAT";
kk = LCP(:,:,i) == 100;                 % indices of all water pixels
landmask(kk) = false;
clear kk i
"Imported land mask from MODIS" %#ok<NOPTS>


% Import current-day snow cover data
% **********************************
%   1)  fix nodata issues (constant night) above arctic circle, set to 100% snowcover (October to
%       March). we ignore the similar problem in Antarctica(!)
%   2)  Set all ocean points to NaN
%   3)  For NaN land points, fill with temporal or spatial nearest neighbors. These missing points
%       are usually over 50% water, which is why they were ignored in the MODIS dataset.
for m = 1 : nmonths
    subdatafname = strcat("F:\IntermediaryFiles\MATfiles\SnowCover\SnowCover_2000-2021_",monthnames(m),".mat");
    load(subdatafname,'snowcover');
    temp = sum(snowcover,2,'omitnan');
    if m < 4 || m > 9
        x = find(temp>0,1,'first') +10; % northernmost non-zero value (w buffer) to address taper issue
        snowcover(1:x,:) = 100; 
    end        
%     snowcover(landmask==0) = NaN;
    snowcover_all(:,:,m) = snowcover./100; % fraction of month-pixel snow covered
end

% b. Gap-fill
for m = 1 : nmonths
    snow = snowcover_all(:,:,m);
%     mis = logical((isnan(snow) + landmask) == 2);
    mis = isnan(snow);
    [a,b] = find(mis);
    for px = 1 : length(a)
        i = a(px);
        j = b(px);
        if m == 1, m1 = 12; else, m1 = m-1; end
        if m == 12, m2 = 1; else, m2 = m+1; end
        val = mean(snowcover_all(i,j,[m1,m2]));
        if isnan(val)
            if j == 1 || j == nlon, jj = [nlon,2]; else, jj = [j-1,j+1]; end 
            if i == 1, ii = 2:3; elseif i == nlat, ii = nlat-2:nlat-1; else, ii = i-1:i+1; end 
            val = mean(snowcover_all(ii,jj,m),'all','omitnan');
            if isnan(val), val = 0; end
        end
        snowcover_all(i,j,m) = val;
        MisSnow(i,j,m) = true;
    end
    strcat("Done importing and filling data for snow cover in",...
        monthnames(m)," (",num2str(length(a))," points)")
end


% 3. Import/resample radiation data
% ---------------------------------
% a. Get data
for m = 1 : nmonths
    subdatafname = ...
        strcat(solarradiation,"Visible_plus_NIR_beam_downward_",lower(monthnames(m)),".rst");
    fid = fopen(subdatafname);
    A = fread(fid,'single');
    direct_solar(:,:,m) = permute(reshape(A,nlon1dg,nlat1dg),[2 1]);
    fclose(fid);
    clear A
    subdatafname = ...
        strcat(solarradiation,"Visible_plus_NIR_diffuse_downward_",lower(monthnames(m)),".rst");
    fid = fopen(subdatafname);
    B = fread(fid,'single');
    diffuse_solar(:,:,m) = permute(reshape(B,nlon1dg,nlat1dg),[2 1]);
    fclose(fid);
    clear B
end
frdiff = diffuse_solar ./ (diffuse_solar + direct_solar);

% b. Resample values in the 0.05 x 0.05 grid (MODIS grid)
scaleratio = 1/latlonscale;
for x = 1 : nlat1dg
    for y = 1 : nlon1dg
        xdestinations = (x-1) *scaleratio +1 : x *scaleratio;
        ydestinations = (y-1) *scaleratio +1 : y *scaleratio;
        for m = 1 : nmonths
            solar_direct_hires(xdestinations,ydestinations,m)  = direct_solar(x,y,m);
            solar_diffuse_hires(xdestinations,ydestinations,m) = diffuse_solar(x,y,m);
            frac_diff_hires(xdestinations,ydestinations,m)     = frdiff(x,y,m);
        end
    end
end
"Imported radiation data" %#ok<NOPTS>

clear direct_solar diffuse_solar frdiff m scaleratio x y


% 4. Import/Resample radiative kernels
% ------------------------------------
for kk = 1 : nkernels
    subdatafname = strcat(radiativekernel,kernelfilenames(kk));
    khires = albkern_hires;
    
    % a. Import kernel data
    if kernels(kk) == "CACK1"
        korig = ncread(subdatafname,kernelvarnames(kk));
        lats = ncread(subdatafname,'Latitude');
        lons = ncread(subdatafname,'Longitude');
    else
        korig = permute(ncread(subdatafname,kernelvarnames(kk)),[2 1 3]);
        lats = ncread(subdatafname,'lat');
        lons = ncread(subdatafname,'lon');
    end
    nkolt = length(lats);
    nkoln = length(lons);
    
    % b. Arrange data in a [NS,WE] array starting with (1,1) at (90N,180W)
    if lons(length(lons)) > 180
        west = lons >= 180;
        east = lons < 180;
        korigwrap = cat(2,korig(:,west,:),korig(:,east,:));
        lonswrap = cat(1,lons(west),lons(east));
    else
        korigwrap = korig;
        lonswrap = lons;
    end
    if lats(1) ~= abs(lats(1))
        korigwrap = flipud(korigwrap);
        lats = flipud(lats);
    end
    
    % d. check for missing data, and fill with 0 when dark
    if sum(isnan(korigwrap),'all') > 0
        for mm = 1 : nmonths
            if sum(isnan(korigwrap(:,:,mm)),'all') > 0
                for lo = 1 : nkoln
                    val = korigwrap(:,lo,mm);
                    if isnan(val(1)) && ismember(mm,[1:3,10:12])
                        ll = find(isnan(val),1,'last');
                        if lats(ll) > 66
                            val(isnan(val)) = 0;
                        else
                            error("code here")
                        end
                    elseif isnan(val(nkolt)) && ismember(mm,4:9)
                        ll = find(isnan(val),1,'first');
                        if lats(ll) < -66
                            val(isnan(val)) = 0;
                        else
                            error("code here")
                        end
                    end
                    korigwrap(:,lo,mm) = val;
                end
            end
        end
    end            
        
    % c. Resample values in the 0.05 x 0.05 grid (MODIS grid) - assuming lat/lon are pixel centers
    ini = 1; leftoveri = 0;
    fsti = round((90 - lats(1)) / latlonscale);
    npixi = zeros(nkolt,1);
    npixj = zeros(nkoln + 1,1);
    for i = 1 : nkolt
        if i < nkolt
            stepi = abs(lats(i+1) - lats(i)) / latlonscale;
            ssti = round(stepi / 2);
            leftoveri = leftoveri + stepi - ssti*2;
            if abs(leftoveri) >= 1
                if leftoveri < 0
                    lasti = ini + fsti + ssti - 2;
                    leftoveri = leftoveri + 1;
                else
                    lasti = ini + fsti + ssti;
                    leftoveri = leftoveri - 1;
                end
            else
                lasti = ini + fsti + ssti - 1;
            end
        else
            lasti = nlat;
        end
        inj = 1;
        fstj = round((lonswrap(1) - 180) / latlonscale);
        leftoverj = ((lonswrap(1) - 180) / latlonscale) - fstj;
        if lons(1) == 0, nbint = nkoln + 1; else, nbint = nkoln; end
        for j = 1 : nbint
            if j == nkoln + 1
                kval = korigwrap(i,1,:);
                lastj = nlon;
                sstj = 0;
            else
                kval = korigwrap(i,j,:);
                if j == nkoln
                    el = lonswrap(1);
                elseif lonswrap(j+1) < lonswrap(j)
                    el = 360 + lonswrap(j+1);
                else
                    el = lonswrap(j+1);
                end
                stepj = (el - lonswrap(j)) / latlonscale;
                sstj = round(stepj / 2);
                leftoverj = leftoverj + stepj - sstj*2;
                if abs(leftoverj) >= 1
                    if leftoverj < 0
                        lastj = inj + fstj + sstj - 2;
                        leftoverj = leftoverj + 1;
                    else
                        lastj = inj + fstj + sstj;
                        leftoverj = leftoverj - 1;
                    end
                else
                    lastj = inj + fstj + sstj - 1;
                end
            end
            for m = 1 : nmonths
                khires(ini:lasti,inj:lastj,m) = kval(m);
            end
            if i == 1
                npixj(j) = length(inj:lastj);
            end
            inj = lastj + 1;
            fstj = sstj;
        end
            
        npixi(i) = length(ini:lasti);
        ini = lasti + 1;
        fsti = ssti;       
    end    
    
    % d. Save kernel
    eval(strcat(kernels(kk)," = khires;"));
    strcat("Done importing kernel data from ",kernels(kk))
    
    clear subdatafname khires korig lats lons nkolt nkoln west east korigwrap ...
    lonswrap ini fsti leftoveri npixi npixj 

    
end
clear kk albkern_hires


% 5. Get pixel areas
% ------------------
%    (in a sphere projection, all pixels on same latitude have same area, so I only need to
%     calculate them along the latitude axis)
earthellipsoid = referenceSphere('earth','m');
pix = zeros(nlat,1);
lat1 = 90;
lon1 = -180;
lon2 = lon1 + latlonscale;
for i = 1 : nlat 
    lat2 = lat1 - latlonscale;
        pix(i) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
    lat1 = lat2;
end
globalpixelarea = repmat(pix,[1 nlon]);     


% 6. Divide all maps in smaller blocks (20*20 deg lat/lon) and save in separate files
% -----------------------------------------------------------------------------------
%    (for computer memory management!)

presland = true(nblocks,1);
for rr = 1 : nblocks
    
    % a. Get coordonates
    y = ceil(rr/nblx);              % block index in the y direction
    x = rr - (y-1)*nblx;            % block index in the x direction
    regj = (x-1)*blocksize + 1 : x * blocksize;
    regi = (y-1)*blocksize + 1 : y * blocksize;
    
    % b. Land cover
    blocklandmask = landmask(regi,regj);
    if sum(blocklandmask,'all') == 0, presland(rr) = false; continue; end
    landcoverprop = LCP(regi,regj,:);
    
    % c. Snow cover
    snowcover = snowcover_all(regi,regj,:);
    
    % d. Radiation
    directsolar = solar_direct_hires(regi,regj,:);
    diffusesolar = solar_diffuse_hires(regi,regj,:);
    diffusefraction = frac_diff_hires(regi,regj,:);
    
    % e. Kernels
    cam3 = CAM3(regi,regj,:);
    cam5 = CAM5(regi,regj,:);
    echam6 = ECHAM6(regi,regj,:);
    hadgem2 = HADGEM2(regi,regj,:);
    cack1 = CACK1(regi,regj,:);
    hadgem3 = HADGEM3(regi,regj,:);
    
    % f. Pixel area
    pixelarea = globalpixelarea(regi,regj);
    
    % g. Save in mat files
    subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
    save(subdatafname,'x','y','regi','regj','landcoverprop','blocklandmask',...
        'snowcover','directsolar','diffusesolar','diffusefraction',...
        'cam3','cam5','echam6','hadgem2','cack1','hadgem3','pixelarea');
    
    strcat("Done with saving data into subregion ",num2str(rr))
    clear x y regi regj landcoverprop blocklandmask snowcover directsolar ...
         diffusesolar diffusefraction cam3 cam5 echam6 hadgem2
end

% h. Save global variables
save(parameterfile,'presland','MisSnow','-append');
save(glodatafname,'presland','globalpixelarea','landmask','MisSnow','-append');

clear CAM3 CAM5 ECHAM6 HADGEM2 solar_direct_hires solar_diffuse_hires frac_diff_hires ...
    snowcover_all LCP landmask



%% C. Import and save albedo atlas data
% *************************************

% 1. Import and save original data
% --------------------------------
for aa = 1 : nalb
    for m = 1 : nmonths       
        albedofname = strcat(albedopaths(aa),"hierarchical_",albid(aa,2),"albedo_igbp_0.05.",...
            albid(aa,1),"_shortwave.",monthnames(m),".mean");
        [varargout] = enviread(albedofname,albedoheader);
        albedo = flipud(varargout);
        clear varargout
        
        for rr = 1 : nblocks
            if presland(rr) == true
                subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
                if m == 1
                    albedoatlas = blockalbedo;
                    load(subdatafname,'regi','regj');
                else
                    load(subdatafname,'albedoatlas','regi','regj');
                end
                temp = permute(albedo(regi,regj,:),[1 2 4 3]);
                temp(temp == nodata_albedo) = NaN;
                temp2 = temp .* scale_albedo;
                albedoatlas(:,:,m,:) = temp2;
                clear temp temp2
                if m == nmonths
                    eval(strcat(albsn(aa)," = albedoatlas;"));
                    varname = cellstr(albsn(aa));
                    save(subdatafname,varname{:},'-append')
                else
                    save(subdatafname,'albedoatlas','-append');
                end
            end
            strcat("done with block ",num2str(rr)," and ",albsn(aa)," - ",monthnames(m))
        end
        clear albedo albedofname rr subdatafname albedoatlas
    end
end

clear aa m blockalbedo


% 2. Fill missing values
% ----------------------
%    (some months/biome-types have no sample anywhere on Earth, snow-covered crops in August for
%     example, so we fill those values with the average of the emcompassing months with data)
albedovars = cellstr(albsn);
for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname,albedovars{:});
        
        for aa = 1 : nalb
            eval(strcat("albedo = ",albsn(aa),";"));
            for bb = 1 : nbiomes
                bioalb = albedo(:,:,:,bb);
                na = squeeze(sum(uint8(isnan(bioalb)),1:2));
                if sum(na) > 1
                    i = find(na > 1);
                    if (i(length(i)) - i(1)) == (length(i) - 1)
                        if i(1)-1 == 0, m1 = 12; else, m1 = i(1)-1; end
                        if i(length(i))+1 == 13, m2 = 1; else, m2 = i(length(i))+1; end
                    elseif i(1) == 1 && i(length(i)) == 12
                        k = find(diff(i)>1);
                        if length(k) == 1
                            m1 = i(k) + 1; m2 = i(k+1) - 1;
                        else
                            error("Code albedo fill values when there is more than one missing period")
                        end
                    else
                        error("Code albedo fill values when there is more than one missing period")
                    end
                    map1 = bioalb(:,:,m1); map2 = bioalb(:,:,m2);
                    bioalb(:,:,i) = repmat(mean(cat(3,map1,map2),3),[1,1,length(i)]);
                end
                albedo(:,:,:,bb) = bioalb;
            end
            eval(strcat(albsn(aa)," = albedo;"));
        end
        
        save(subdatafname,albedovars{:},'-append')
    end
end

clear(albedovars{:})
clear rr subdatafname aa albedo bb bioalb na i k m1 m2 map1 map2



