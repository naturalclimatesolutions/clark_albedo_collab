% NCS - Global Albedo - TNC-Clark Collaboration - Glasgow COP results
%
% This code is used to produce data for the Glasgow COP meeting, as a first estimate of albedo
%   offsets due to land cover change. Used data are MODIS/Terra for land-cover mask (we define  
%   values everywhere, independent of plausible presence of the studied land cover),
%   the Gao et al. 2014 albedo catalog, the present-day snow cover (from MODIS), and the CAM3 
%   (Shell et al. 2008) CAM5, ECHAM6 and HADGEM2 radiative kernels. Those are all "hard-coded".
%
% Note: Kernels do not specify if lat/lon coordinate correspond to the center of the pixel or the
%       western/southern boundary. I assumed they are the center points.
% Note: For any land pixel in MODIS/Terra (at least 1% non-water), we calculate albedo impact of any
%       land cover conversion regardless of plausibility. Even though I initially calculated those
%       values for Antarctica as well, I eliminated all values South of 60S for readability of the
%       output maps.
%
% Created 6/30/2021 by Natalia Hasler, Clark University

% Last modified: NH 9/27/2021
%                (See Modification history at end of file)



clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A.  Define inputs/outputs
% **************************

% **************************************
% 1. User-defined variables
% -------------------------
computer = "GEO-005199";                % where the routine is run (for directory paths):
                                        %   GEO-007264 (my laptop)
                                        %   GEO-005199 (old supercomputer)
                                        %   GEO-006381 (new supercomputer)

pathways = ...                          % name of studied pathways
    ["deforestation to croplands";"deforestation to grasslands";...
    "deforestation to urban developement";"reforestation from open shrublands";...
    "reforestation from closed shrublands";"grassland to open savanna";"grassland to woody savanna"];
pwnames = ["Forest2Crop";...            % pathway names to use as folder names
    "Forest2Grass";"Forest2Urban";...
    "OpenShrub2Forest";"ClosedShrub2Forest";"Grass2Savanna";"Grass2WoodySavanna"];
npw = length(pathways);
forestbiomes = ["ENF","EBF","DNF","DBF","MF"];
nfrst = length(forestbiomes);

plc = strings(nfrst,3,npw-2);           % String array of 3-letter acronyms of old/new land cover in
                                        %   the desired pathway (one/line) with 3 columns as
                                        %   1.  3-letter acronym of initial land cover
                                        %   2.  3-letter acronym of final land cover
                                        %   3.  pathway number (matching line number in "pathways")
plc(:,1,1:3) = repmat(forestbiomes,[1,1,3]);
plc(:,2,1) = repmat("CRO",[nfrst,1]);
plc(:,2,2) = repmat("GRA",[nfrst,1]);
plc(:,2,3) = repmat("URB",[nfrst,1]);
plc(:,2,4:5) = repmat(forestbiomes',[1,1,2]);
plc(:,1,4) = repmat("OSH",[nfrst,1]);
plc(:,1,5) = repmat("CSH",[nfrst,1]);
plc(:,3,:) = permute(repmat(1:5,[nfrst,1]),[1,3,2]);
pathwayslandcover = reshape(permute(plc,[1,3,2]),[nfrst*(npw-2),3]);
pathwayslandcover = cat(1,pathwayslandcover,["GRA","SAV","6"],["GRA","WSA","7"]);

% ***************************************


% 2. Define file/folder locations/names
% -------------------------------------
switch computer                         % directory paths for project data (depending on computer)
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
        mapdir = strcat('G:\GlobalAlbedo\');
        resultdir = rootdir;
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
        mapdir = strcat('G:\TNC\GlobalAlbedo\');
        resultdir = 'E:\NCS-GlobalAlbedo\';
end

% b. Input dataset files/folder paths
landcovermap = strcat(mapdir,"LandCover\MCD12C1.A2016001.006.2018324172410.hdf");
albedosnowfree = strcat(mapdir,"AlbedoAtlas\snow_free_hierarchical_v3\");
albedosnowcovered = strcat(mapdir,"AlbedoAtlas\snow_covered_hierarchical_v3\");
albedoheader = strcat(albedosnowcovered,"hierarchical_snow_albedo_igbp_0.05.bsa_shortwave.Jan.hdr");
albedopaths = cat(1,albedosnowfree,albedosnowfree,albedosnowcovered,albedosnowcovered);
snowcoverdata = strcat(mapdir,"SnowData\");
solarradiation = strcat(mapdir,"solarflux\");
radiativekernel = strcat(mapdir,"RadiativeKernels\");
regoutputfiles = strcat(rootdir,"SubData\");
boxpath = "C:\Users\nhasler\Box\Albedo (Clark-TNC)\";
resultslocalfolder = strcat(resultdir,"GlasgowCOP\");


% 3. Define input parameters
% --------------------------
Ca_PgC  = 400*2.13;                     % current base CO2 in atmos [Pg C]
Aglobe  = 5.1007e14;                    % global surface area [m2]
latlonscale = 0.05;                     % MODIS GCM scale

IGBPBiomes = ["Water bodies"; ...       % IGBP Biomes names
    "Evergreen needleleaf forests"; "Evergreen broadleaf forests"; "Deciduous needleleaf forests";...
    "Deciduous broadleaf forests"; "Mixed forests"; "Closed shrublands"; "Open shrublands"; ...
    "Woody savannas"; "Savannas"; "Grasslands"; "Permanent wetlands"; "Croplands"; ...
    "Urban and built-up lands"; "Cropland/natural vegetation mosaics"; "Snow and ice"; "Barren"];
biomes_igbp = ["WAT","ENF","EBF",...    % IGBP Biomes 3-letter acronyms
    "DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","MOS","SNO","BAR"];
nbiomes = length(biomes_igbp);          % number of IGBP biome categories

monthnames = ["Jan","Feb","Mar",...     % months 3-letter names (used in some file names)
    "Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
nmonths = length(monthnames);           % number of months

albedo_types = ...                      % Albedo values found in albedo atlas
    ["snow-free black-sky","snow-free white-sky","snow-covered black-sky","snow-covered white-sky"];
albsn = ["nosnowblacksky",...           % Albedo short names (used as variable names)
    "nosnowwhitesky","snowcovblacksky","snowcovwhitesky"];
albid = cat(2,...                       % String values used in albedo file naming
    ["bsa";"wsa";"bsa";"wsa"],["";"";"snow_";"snow_"]);
nodata_albedo = 2^15-1;                 % Albedo "nodata" value
scale_albedo = 0.001;                   % Albedo scale value
nalb = length(albedo_types);            % number of albedo types

kernels = ["CAM3","CAM5",...            % Radiative kernel names
    "ECHAM6","HADGEM2","CACK1","HADGEM3"];
nkernels = length(kernels);             % number of radiative kernels

nlcc = length(pathwayslandcover);       % number of individual land cover conversions
pwid = ...                              % correspondance between individual land cover conversions 
    str2double(pathwayslandcover(:,3)); %   and pathway groups

nlat1dg = 180;                          % 1 degree latitudes
nlon1dg = 360;                          % 1 degree longitudes
nlat = nlat1dg ./ latlonscale;          % number of pixels in MODIS-GCM latitude-wise
nlon = nlon1dg ./ latlonscale;          % number of pixels in MODIS-GCM longitude-wise

kernelfilenames = strings(nkernels,1);  % kernel file names array (used in read loop)
kernelvarnames = strings(nkernels,1);   % kernel variable name in original file
kernelscale = ones(nkernels,1);         % values of radiaive forcings are often stored for 0.01
                                        % albedo changes, and sometimes sign is opposite too ...
for kk = 1 : nkernels
    kname = kernels(kk);
    switch kname
        case "CAM3"
            kernelfilenames(kk) = "CAM3\CAM3_albedo_sw_kernel.nc";
            kernelvarnames(kk) = "monkernel";
            kernelscale(kk) = 100;
        case "CAM5"
            kernelfilenames(kk) = "CAM5\alb.kernel.nc";
            kernelvarnames(kk) = "FSNT";
            kernelscale(kk) = 100;
        case "ECHAM6"
            kernelfilenames(kk) = "ECHAM6\ECHAM6_CTRL_kernel.nc";
            kernelvarnames(kk) = "A_srad0";
            kernelscale(kk) = 100;
        case "HADGEM2"
            kernelfilenames(kk) = "HADGEM2\HadGEM2_sw_TOA_L38.nc";
            kernelvarnames(kk) = "albedo_sw";
            kernelscale(kk) = 100;
        case "CACK1"
            kernelfilenames(kk) = "CACKv1.0\CACKv1.0.nc";
            kernelvarnames(kk) = "CACK CM";
            kernelscale(kk) = -1;
        case "HADGEM3"
            kernelfilenames(kk) = "HadGEM3-GA7.1_TOA_kernel_L85.nc";
            kernelvarnames(kk) = "albedo_sw";
            kernelscale(kk) = 100;
        otherwise
            error("This kernel as not been defined yet!")
    end
end


% 4. Define storage arrays
% ------------------------
% a. Global input arrays
snowcover_all = zeros(nlat,nlon,nmonths);
direct_solar = zeros(nlat1dg,nlon1dg,nmonths);
diffuse_solar = zeros(nlat1dg,nlon1dg,nmonths);
solar_direct_hires  = zeros(nlat,nlon,nmonths);
solar_diffuse_hires = zeros(nlat,nlon,nmonths);
frac_diff_hires     = zeros(nlat,nlon,nmonths);
albkern_hires = zeros(nlat,nlon,nmonths);

% b. Regional input arrays
%    (for memory reasons, I have to split the planet in 20X20 degrees lat/lon chunks)
blocksize = 20 / latlonscale;
nbly = nlat / blocksize;
nblx = nlon / blocksize;
nblocks = nblx * nbly;
blockalbedo = zeros(blocksize,blocksize,nmonths,nbiomes);

% c. Global output arrays
forestpathways = NaN(nlat,nlon,nkernels,nfrst); %#ok<NASGU>
nonforestpathways = NaN(nlat,nlon,nkernels); %#ok<NASGU>

for pp = 1 : npw
    if sum(pwid==pp) == nfrst
        eval(strcat(pwnames(pp)," = forestpathways;"))
        eval(strcat("Miss",pwnames(pp)," = missfpw;"))
    elseif sum(pwid==pp) == 1
        eval(strcat(pwnames(pp)," = nonforestpathways;"))
        eval(strcat("Miss",pwnames(pp)," = missnfpw;"))
    else
        error(strcat("The ",pathways(pp)," pathway has more than one land cover change, ",...
            "please modify code in section A.3.c"))
    end
end
        
MisSnow = false(nlat,nlon,nmonths);     % Used to record where I gap-filled the snow data

clear forestpathways nonforestpathways missnfpw missfpw

% d. Regional output arrays
regfCO2e= NaN(blocksize,blocksize,nkernels,nfrst);
regnfCO2e= NaN(blocksize,blocksize,nkernels);
regfmiss = false(blocksize,blocksize,nmonths,nfrst);
regnfmiss = false(blocksize,blocksize,nmonths);


% 5. Save input variables for later use
% -------------------------------------
glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
save(glodatafname,'regoutputfiles','boxpath','nblocks','resultslocalfolder','blocksize',...
    'nblx','nbly','nlon','nlat','nkernels','IGBPBiomes','biomes_igbp','Ca_PgC',...
    'nfrst','monthnames','albedo_types','albsn','albid','fillvalue_albedo','scale_albedo','nalb',...
    'kernels','kernelscale','nbiomes','nmonths','Aglobe','latlonscale',...
    'pathways','pwnames','npw','pathwayslandcover','nlcc','pwid');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B. Import, fill/correct and resample data (except albedo)
% ******************************************

% 1. Import land cover data
% -------------------------
% Note: IGBP classification is given as "land cover type 1"
%       I import the land cover proportion for each land cover-type ("Land_Cover_Type_1_Percent") to
%       calculate albedo-induced radiative forcing on every land pixel that has some land cover.
% Note: MODIS/Terra grid is in 0.05 x 0.05 degree grid, with upperleft corner (90,-180) lat/lon
% a. Import data
LCP = hdfread(landcovermap,"Land_Cover_Type_1_Percent");
% b. Find land mask (any pixel that has at least 1% non-water cover)
landmask = true(nlat,nlon);
i = biomes_igbp =="WAT";
kk = LCP(:,:,i) == 100;                 % indices of all water pixels
landmask(kk) = false;
clear kk i
"Imported land mask from MODIS" %#ok<NOPTS>


% 2. Import current-day snow cover data
% -------------------------------------
% Note: Snow Cover climatology has been pre-calculated from MODIS/Terra Snow Cover Monthly L3
%       Global 0.05Deg CMG, Version 6 with script "generate_snow_cover_climatology.m"
% Additional issues fixed here:
%   1)  fix nodata issues (constant night) above arctic circle, set to 100% snowcover (October to
%       March). we ignore the similar problem in Antarctica(!)
%   2)  Set all ocean points to NaN
%   3)  For NaN land points, fill with temporal or spatial nearest neighbors. These missing points
%       are usually over 50% water, which is why they were ignored in the MODIS dataset.

% a. Import data
for m = 1 : nmonths
    subdatafname = strcat(snowcoverdata,"SnowCover_2000-2021_",monthnames(m),".mat");
    load(subdatafname,'snowcover');
    temp = sum(snowcover,2,'omitnan');
    if m < 4 || m > 9
        x = find(temp>0,1,'first') +10; % northernmost non-zero value (w buffer) to address taper issue
        snowcover(1:x,:) = 100; 
    end        
    snowcover(landmask==0) = NaN;
    snowcover_all(:,:,m) = snowcover./100; % fraction of month-pixel snow covered
end

% b. Gap-fill
for m = 1 : nmonths
    snow = snowcover_all(:,:,m);
    mis = logical((isnan(snow) + landmask) == 2);
    [a,b] = find(mis);
    for px = 1 : length(a)
        i = a(px);
        j = b(px);
        if m == 1, m1 = 12; else, m1 = m-1; end
        if m == 12, m2 = 1; else, m2 = m+1; end
        val = mean(snowcover_all(i,j,[m1,m2]));
        if isnan(val)
            if j == 1 || j == nlon, jj = [nlon,2]; else, jj = [j-1,j+1]; end %#ok<BDSCI>
            if i == 1, ii = 2:3; elseif i == nlat, ii = nlat-2:nlat-1; else, ii = i-1:i+1; end %#ok<BDSCI>
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
    end
    if lats(1) ~= abs(lats(1))
        korigwrap = flipud(korigwrap);
        lats = flipud(lats);
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
    
end
clear kk subdatafname khires albkern_hires korig lats lons nkolt nkoln west east korigwrap ...
    lonswrap ini fsti leftoveri npixi npixj 


% 5. Get pixel areas
% ------------------
%    (in a sphere projection, all pixels on same latitude have same area, so I only need to
%     calculate them along the latitude axis)
earthellipsoid = referenceSphere('earth','m');
pix = zeros(nlat,1);
lat1 = 90;
lon1 = -180;
lon2 = lon1 + latlonscale;
for i = 1 : nlat %#ok<BDSCI>
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
save(glodatafname,'presland','globalpixelarea','landmask','MisSnow','-append');

clear CAM3 CAM5 ECHAM6 HADGEM2 solar_direct_hires solar_diffuse_hires frac_diff_hires ...
    snowcover_all LCP landmask




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% D. Compute TOA RF and CO2e of land conversion scenarios
% ********************************************************

% 1. Loop through sub-regions
% ---------------------------
tilediv = 4;
tlsize = blocksize / tilediv;
nbtiles = tilediv * tilediv;
        
for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname);
        
        
        % 2. Sub-divide region in smaller tiles for use in parfor loop
        % ------------------------------------------------------------
        tilemvars = NaN(tlsize,tlsize,nmonths,nbtiles);
        tilealvar = NaN(tlsize,tlsize,nmonths,nbiomes,nbtiles);
        tileco2e = NaN(tlsize,tlsize,nkernels,nlcc,nbtiles);
        tilemiss = false(tlsize,tlsize,nmonths,nlcc,nbtiles);
        
        tlnsbs = tilealvar; tlnsws = tilealvar; tlscbs = tilealvar; tlscws = tilealvar;
        tlsnowcover = tilemvars; tldiffrac = tilemvars;
        
        tlkernels = NaN(tlsize,tlsize,nmonths,nkernels,nbtiles);
        tllandmask = false(tlsize,tlsize,nbtiles);
        tlpixarea = NaN(tlsize,tlsize,nbtiles);
                
        clear tilemvars tilealvar
        
        for tl = 1 : nbtiles
            y = ceil(tl/tilediv);
            x = tl - (y-1)*tilediv;
            tlj = (x-1)*tlsize + 1 : x * tlsize;
            tli = (y-1)*tlsize + 1 : y * tlsize;
            
            tlnsbs(:,:,:,:,tl) = nosnowblacksky(tli,tlj,:,:);
            tlnsws(:,:,:,:,tl) = nosnowwhitesky(tli,tlj,:,:);
            tlscbs(:,:,:,:,tl) = snowcovblacksky(tli,tlj,:,:);
            tlscws(:,:,:,:,tl) = snowcovwhitesky(tli,tlj,:,:);
            
            tlsnowcover(:,:,:,tl) = snowcover(tli,tlj,:);
            tldiffrac(:,:,:,tl) = diffusefraction(tli,tlj,:);
            
            tllandmask(:,:,tl) = blocklandmask(tli,tlj);
            tlpixarea(:,:,tl) = pixelarea(tli,tlj);
            
            for kk = 1 : nkernels
                eval(strcat("kern = ",lower(kernels(kk)),";"))
                tlkernels(:,:,:,kk,tl) = kern(tli,tlj,:);
            end
            
            clear x y tli tlj kk
        end

        
        % 3. Loop through tiles, define arrays and find land pixels
        % ---------------------------------------------------------
        parfor tl = 1 : nbtiles
            
            sm_nosnow_blacksky = tlnsbs(:,:,:,:,tl);
            sm_nosnow_whitesky = tlnsws(:,:,:,:,tl);
            sm_snowcov_blacksky = tlscbs(:,:,:,:,tl);
            sm_snowcov_whitesky = tlscws(:,:,:,:,tl);
            
            sm_snowcover = tlsnowcover(:,:,:,tl);
            sm_diffuse_fraction = tldiffrac(:,:,:,tl);
            sm_kernels = tlkernels(:,:,:,:,tl);
            
            sm_landmask = tllandmask(:,:,tl);
            sm_pixelarea = tlpixarea(:,:,tl);
            
            sm_pathwayco2 = NaN(tlsize,tlsize,nkernels,nlcc);
            sm_pathwaymiss = false(tlsize,tlsize,nmonths,nlcc);
            
            [a,b] = find(sm_landmask);
        
        
            % 4. Loop through land pixels/pathways and calculate monthly radiative forcings
            % -----------------------------------------------------------------------------
            for pt = 1 : length(a)
                i = a(pt);
                j = b(pt);
                
                % a. Loop through individual pathways
                for ll = 1 : nlcc
                    lc1 = biomes_igbp == pathwayslandcover(ll,1); %#ok<PFBNS>
                    lc2 = biomes_igbp == pathwayslandcover(ll,2);
                    ta = zeros(nmonths,nkernels);
                    
                    % b. Loop through months                               
                    for m = 1 : nmonths
                    
                        % c. Calculate surface albedo differences
                        daNSBS = sm_nosnow_blacksky(i,j,m,lc1) - sm_nosnow_blacksky(i,j,m,lc2);
                        daNSWS = sm_nosnow_whitesky(i,j,m,lc1) - sm_nosnow_whitesky(i,j,m,lc2);
                        daSCBS = sm_snowcov_blacksky(i,j,m,lc1) - sm_snowcov_blacksky(i,j,m,lc2);
                        daSCWS = sm_snowcov_whitesky(i,j,m,lc1) - sm_snowcov_whitesky(i,j,m,lc2);

                        % d. Calculate monthly Surface Weighted Albedo Differences
                        %    (to avoid NaN values due to undefined albedos in irrelevant situations,
                        %     I add a cumbersome if statement ... now that we filled these, it
                        %     should not be necessary, but I left it as is)
                        snow = sm_snowcover(i,j,m);
                        if snow == 0
                            Wda = squeeze(daNSBS .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                                daNSWS .* sm_diffuse_fraction(i,j,m));
                        elseif snow == 1
                            Wda = squeeze(daSCBS .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                                daSCWS .* sm_diffuse_fraction(i,j,m));
                        else
                            Wda = squeeze(...
                                daNSBS .* (1 - snow) .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                                daNSWS .* (1 - snow) .* sm_diffuse_fraction(i,j,m) + ...
                                daSCBS .* snow .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                                daSCWS .* snow .* sm_diffuse_fraction(i,j,m));
                        end
                                        
                        % e. Compute monthly TOA RF in [W m-2]
                        for kk = 1 : nkernels
                            ta(m,kk) =  Wda .* sm_kernels(i,j,m,kk) .*kernelscale(kk); %#ok<PFBNS>
                        end
                    
                    end     % end of months loop
                
                
                    % 5. Save missing values for quality control/improvements
                    % -------------------------------------------------------
                    %    (Kernels don't have missing values, so I can ignore that dimension)
                    sm_pathwaymiss(i,j,:,ll) = isnan(squeeze(ta(:,1)));
                
                
                    % 6. Calculate CO2 equivalent
                    % ---------------------------
                    % a. Annual mean RF
                    arf = squeeze(mean(ta,1,'omitnan'));
                    
                    % b. Global annual mean RF
                    glob = arf .* sm_pixelarea(i,j) ./ Aglobe;
                
                    % c. CO2e (Pg C)
                    %    (also, flip sign for ease of read)
                    PgC_TOA = Ca_PgC .* exp(-1 .* glob ./ 5.35);
                                    
                    % d. in tonnes CO2 ha-2
                    CO2eq_TOA = (PgC_TOA - Ca_PgC) .* 1e13 ./ sm_pixelarea(i,j) .* 44/12;
                
                    % e. Save values in tile arrays
                    sm_pathwayco2(i,j,:,ll) = CO2eq_TOA;
                
                end     % end loop on all individual pathways
                                
            end     % end loop on all forested pixels
            
            
            % 8. Save tile results back into regional arrays
            % ----------------------------------------------
            tileco2e(:,:,:,:,tl) = sm_pathwayco2;
            tilemiss(:,:,:,:,tl) = sm_pathwaymiss;
            
            strcat("done with RF/CO2e calculation for subregion ",num2str(rr)," - tile #",...
                num2str(tl)," (",num2str(length(a))," points)")

        end     % end parfor on tiles
        
        
        % 9. Recombine tiles into one regional arrays
        % -------------------------------------------        
        for pp = 1 : npw
            if sum(pwid==pp) == nfrst
                eval(strcat("regco2",lower(pwnames(pp))," = regfCO2e;"));
                eval(strcat("regmiss",lower(pwnames(pp))," = regfmiss;"));
            else
                eval(strcat("regco2",lower(pwnames(pp))," = regnfCO2e;"));
                eval(strcat("regmiss",lower(pwnames(pp))," = regnfmiss;"));
            end
        end
        
        if ismember(rr,canblocks)
            canrf = regcanmonths; canda = regcanmonths; cangra = regcan; cancro = regcan;
        end

        for tl = 1 : nbtiles
            y = ceil(tl/tilediv);
            x = tl - (y-1)*tilediv;
            tlj = (x-1)*tlsize + 1 : x * tlsize;
            tli = (y-1)*tlsize + 1 : y * tlsize;
            
            for pp = 1 : npw
                ii = pwid == pp;
                eval(strcat("regco2",lower(pwnames(pp)),"(tli,tlj,:,:) = tileco2e(:,:,:,ii,tl);"));
                eval(strcat("regmiss",lower(pwnames(pp)),"(tli,tlj,:,:) = tilemiss(:,:,:,ii,tl);"));
            end
        end            
        
        
        % 10. Save values in global arrays
        % --------------------------------
        for pp = 1 : npw
            eval(strcat(pwnames(pp),"(regi,regj,:,:) = regco2",lower(pwnames(pp)),";"));
            eval(strcat("Miss",pwnames(pp),"(regi,regj,:,:) = regmiss",lower(pwnames(pp)),";"));
        end

        regvars = cellstr(cat(1,strcat("regco2",lower(pwnames)),...
            strcat("regmiss",lower(pwnames))));
        save(subdatafname,regvars{:},'-append')
        
        strcat("done with RF/CO2e calculation for subregion ",num2str(rr),...
            " (",num2str(sum(blocklandmask,'all'))," points)")
        
        clear(regvars{:})
        clear nosnowblacksky nosnowwhitesky snowcovblacksky snowcovwhitesky ...
            x y regi regj landcoverprop blocklandmask snowcover directsolar ...
            diffusesolar diffusefraction cam3 cam5 echam6 hadgem2 ...
            canrf canda cangra cancro tlnsbs tlnsws tlscbs tlscws tlsnowcover tldiffrac ...
            tlkernels tllandmask tlpixarea tileco2e tilemiss ...
            tlcanrf tlcanad tlcangra tlcancro

        
    end     % end if region has land data
end
        
datafname = strcat(regoutputfiles,"AlbedoFinalData6ker.mat");
globalvars = cellstr(cat(1,pwnames,strcat("Miss",pwnames)));
save(datafname,globalvars{:},'-v7.3')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% E. Export arrays in separate geotiff maps
% ******************************************
%    (I also eliminate Antarctica here, should have done it earlier!)

% 1. Define Antarctica (or more precisely South of 60S)
% --------------------
lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[~,lats] = meshgrid(lon05,lat05);
latA = -60;
Antarctica = lats <= latA;


% 2. Save data in GeoTIFF fornat
% ------------------------------

bbox = [-180,-90;180,90];

for pp = 1 : npw
    pname = pwnames(pp);
    ii = find(pwid==pp);
    for ll = 1 : length(ii)
        z = ii(ll);
        for kk = 1 : nkernels
            eval(strcat("data = ",pname,"(:,:,kk,ll);"))
            data(Antarctica) = NaN; %#ok<SAGROW>
            data(landmask==0) = NaN; %#ok<SAGROW> % this should not be necessary, but somehow I still got zero values
            fname = strcat(resultslocalfolder,pname,"\",pathwayslandcover(z,1),"2",...
                pathwayslandcover(z,2),"_",kernels(kk),".tif");
            geotiffwrite_easy(fname,bbox,data);
        end
    end
    clear pname ii ll z data fname
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MODIFICATIONS
% **************
% (i don't list any cosmetic changes)

% 9/27/21 NH Added the CACKv1.0 and HadGEM3 kernels
% 9/10/21 NH Gap-filled albedo data for missing months
% 9/8/21 NH Made pathways more general to be able to add/remove some without having to change the
%           code.
% 9/2/21 NH Modified snow read
% 9/2/21 NH Eliminated NaN due to missing snow albedo for month/pixels with no snow
% 8/30/21 NH Output as geotiff maps
% 8/26/21 NH Changed calculation to include all land pixels and added parfor loop
% 8/26/21 NH Changed units to tCO2/ha
    
