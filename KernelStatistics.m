% Global Albedo - Kernel Statistics
%
% Uses results from the GlobalAlbedo Glasgow COP runs, namely outputs from different kernels data
%   and compiles some statistics. It also computes statistics on the original input kernels data.
% Note: To avoid memory issues, results are kept in submaps
%
% Created 10/19/2021 by Natalia Hasler, Clark University


clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Define and load variables
% *****************************


% **********************
% 1. User-Defines variables
% -------------------------

computer = "GEO-005199";                % where the routine is run (for directory paths):
                                        %   GEO-007264 (my laptop)
                                        %   GEO-005199 (old supercomputer)
                                        %   GEO-006381 (new supercomputer)
desiredstatistics = ...                 % list of statistics to be performed
    ["min","max","mean","median","range","standard deviation"];
thresholds = [0.1,0.25,0.5,1];          % Inclusion threshold: see where each kernels falls within
                                        %   ensemble mean +/- thresholds

% **********************


% 2. Define inputs
% ----------------
switch computer
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
end

glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'IGBPBiomes','biomes_igbp','latlonscale','nblocks','nlat','nlon','nlcc',...
    'nmonths','nkernels','npw','presland','pathwayslandcover','pwnames','pwid','kernels',...
    'regoutputfiles','resultslocalfolder','landmask','blocksize','kernelscale')

monthln = ["January","February","March","April","May","June","July","August",...
    "September","October","November","December"];

regvars = cellstr(cat(1,strcat("regco2",lower(pwnames))));
kervars = cellstr(lower(kernels));
nstats = length(desiredstatistics);
nt = length(thresholds);

statisticsmaps = NaN(nlat,nlon,nstats);

expectedstat = ["min","max","mean","median","range","standard deviation"];
for ss = 1 : nstats
    stat = desiredstatistics(ss);
    if ismember(stat,["min","max","mean","median","range","standard deviation"]) == false
        error(strcat("Your desired statistic '",stat,"' does not belong to the expected list: ",...
            join(expectedstat,", "),". Check spelling or add code if new statistic is desired!"));
    end
end
clear expectedstat ss stat




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B. Compute statistics on different kernel outcomes
% ***************************************************

% 1. Loop through sub-regions (and get/rearrange input data)
% ---------------------------
tilediv = 4;
tlsize = blocksize / tilediv;
nbtiles = tilediv * tilediv;
        
for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname,regvars{:},kervars{:},'blocklandmask');
        
        smallkernels = NaN(blocksize,blocksize,nmonths,nkernels);
        regco2 = NaN(blocksize,blocksize,nkernels,nlcc);
        for kk = 1 : nkernels
        	eval(strcat("smallkernels(:,:,:,kk) = ",lower(kernels(kk))," .* kernelscale(kk);"));
        end
        for pp = 1 : npw
            ii = pwid == pp; %#ok<NASGU>
            eval(strcat("regco2(:,:,:,ii) = regco2",lower(pwnames(pp)),";"));
            
        end
        
        regstats = NaN(blocksize,blocksize,nstats,nlcc);
        kmstat = NaN(blocksize,blocksize,nmonths,nstats);
        kanmestat = NaN(blocksize,blocksize,nstats);
        kalmostat = NaN(blocksize,blocksize,nstats);
        kerdev = false(blocksize,blocksize,nmonths,nt,nkernels);
        kerinstd = false(blocksize,blocksize,nmonths,nkernels);

        
        % 2. Sub-divide region in smaller tiles for use in parfor loop
        % ------------------------------------------------------------
        tileco2e = NaN(tlsize,tlsize,nkernels,nlcc,nbtiles);
        tilekernels = NaN(tlsize,tlsize,nmonths,nkernels,nbtiles);
        tileco2stat = NaN(tlsize,tlsize,nstats,nlcc,nbtiles);
        tilekmstat = NaN(tlsize,tlsize,nmonths,nstats,nbtiles);
        tilekeranmestat = NaN(tlsize,tlsize,nstats,nbtiles);
        tilekeralmostat = NaN(tlsize,tlsize,nstats,nbtiles);
        tilekerdev = false(tlsize,tlsize,nmonths,nt,nkernels,nbtiles);
        tilekerinstd = false(tlsize,tlsize,nmonths,nkernels,nbtiles);
        tilelandmask = false(tlsize,tlsize,nbtiles);
                
        for tl = 1 : nbtiles
            y = ceil(tl/tilediv);
            x = tl - (y-1)*tilediv;
            tlj = (x-1)*tlsize + 1 : x * tlsize;
            tli = (y-1)*tlsize + 1 : y * tlsize;
            
            tileco2e(:,:,:,:,tl) = regco2(tli,tlj,:,:);
            tilekernels(:,:,:,:,tl) = smallkernels(tli,tlj,:,:);
            tilelandmask(:,:,tl) = blocklandmask(tli,tlj);
            
            clear x y tli tlj
        end

        
        % 3. Loop through tiles, define arrays and find land pixels
        % ---------------------------------------------------------
        parfor tl = 1 : nbtiles
            
            sm_co2 = tileco2e(:,:,:,:,tl);
            sm_kernels = tilekernels(:,:,:,:,tl);
            sm_co2stat = tileco2stat(:,:,:,:,tl);
            sm_kermostat = tilekmstat(:,:,:,:,tl);
            sm_keranmeanstat = tilekeranmestat(:,:,:,tl);
            sm_kerallmonthstat = tilekeralmostat(:,:,:,tl);
            sm_kerdev = tilekerdev(:,:,:,:,:,tl);
            sm_kerinstd = tilekerinstd(:,:,:,:,tl);
            sm_landmask = tilelandmask(:,:,tl);
            
            [a,b] = find(sm_landmask);
            
        
            % 4. Loop through land pixels and calculate statistics
            % ----------------------------------------------------
            for pt = 1 : length(a)
                i = a(pt);
                j = b(pt);
                
                co2 = NaN(nstats,nlcc);
                kermonth = NaN(nmonths,nstats);
                kerannualmean = NaN(1,nstats);
                kerallmonth = NaN(1,nstats);
                
                kernannualmean = mean(sm_kernels(i,j,:,:),3,'omitnan');
                
                for ss = 1 : nstats
                    stat = desiredstatistics(ss); %#ok<PFBNS>
                    switch stat
                        case "min"
                            co2(ss,:) = min(sm_co2(i,j,:,:),[],3,'omitnan');
                            kermonth(:,ss) = min(sm_kernels(i,j,:,:),[],4,'omitnan');
                            kerannualmean(ss) = min(kernannualmean,[],4,'omitnan');
                            kerallmonth(ss) = min(sm_kernels(i,j,:,:),[],'all','omitnan');
                        case "max"
                            co2(ss,:) = max(sm_co2(i,j,:,:),[],3,'omitnan');
                            kermonth(:,ss) = max(sm_kernels(i,j,:,:),[],4,'omitnan');
                            kerannualmean(ss) = max(kernannualmean,[],4,'omitnan');
                            kerallmonth(ss) = max(sm_kernels(i,j,:,:),[],'all','omitnan');                                                   
                        case "median"
                            co2(ss,:) = median(sm_co2(i,j,:,:),3,'omitnan');
                            kermonth(:,ss) = median(sm_kernels(i,j,:,:),4,'omitnan');
                            kerannualmean(ss) = median(kernannualmean,4,'omitnan');
                            kerallmonth(ss) = median(sm_kernels(i,j,:,:),'all','omitnan');
                        case "mean"
                            co2(ss,:) = mean(sm_co2(i,j,:,:),3,'omitnan');
                            kermonth(:,ss) = mean(sm_kernels(i,j,:,:),4,'omitnan');
                            kerannualmean(ss) = mean(kernannualmean,4,'omitnan');
                            kerallmonth(ss) = mean(sm_kernels(i,j,:,:),'all','omitnan');
                        case "standard deviation"
                            co2(ss,:) = std(sm_co2(i,j,:,:),0,3,'omitnan');
                            kermonth(:,ss) = std(sm_kernels(i,j,:,:),0,4,'omitnan');
                            kerannualmean(ss) = std(kernannualmean,0,4,'omitnan');
                            kerallmonth(ss) = std(sm_kernels(i,j,:,:),0,'all','omitnan');
                        case "range"
                            co2(ss,:) = range(sm_co2(i,j,:,:),3);
                            kermonth(:,ss) = range(sm_kernels(i,j,:,:),4);
                            kerannualmean(ss) = range(kernannualmean,4);
                            kerallmonth(ss) = range(sm_kernels(i,j,:,:),'all');                           
                    end
                end
                
                sm_co2stat(i,j,:,:) = co2;
                sm_kermostat(i,j,:,:) = kermonth;
                sm_keranmeanstat(i,j,:) = kerannualmean;
                sm_kerallmonthstat(i,j,:) = kerallmonth;
                
                imean = desiredstatistics == "mean";
                istd = desiredstatistics == "standard deviation";
                for kk = 1 : nkernels
                    dev = abs(squeeze(sm_kernels(i,j,:,kk)) - kermonth(:,imean));
                    
                    for m = 1 : nmonths
                        ii = find(thresholds*100 > dev(m),1,'first');
                        sm_kerdev(i,j,m,ii,kk) = true;
                        if dev(m) <= kermonth(m,istd)
                            sm_kerinstd(i,j,m,kk) = true;
                        end
                    end
                end
                
            end     % end loop on land pixels
            
            
            % 5. Save tile results back into regional arrays
            % ----------------------------------------------
            tileco2stat(:,:,:,:,tl) = sm_co2stat;
            tilekmstat(:,:,:,:,tl) = sm_kermostat;
            tilekeranmestat(:,:,:,tl) = sm_keranmeanstat;
            tilekeralmostat(:,:,:,tl) = sm_kerallmonthstat;
            tilekerdev(:,:,:,:,:,tl) = sm_kerdev;
            tilekerinstd(:,:,:,:,tl) = sm_kerinstd;
            
            strcat("done statistics calculation for subregion ",num2str(rr)," - tile #",...
                num2str(tl)," (",num2str(length(a))," points)")

        end     % end parfor on tiles
        
        
        % 6. Recombine tiles into one regional arrays
        % -------------------------------------------        
        for tl = 1 : nbtiles
            y = ceil(tl/tilediv);
            x = tl - (y-1)*tilediv;
            tlj = (x-1)*tlsize + 1 : x * tlsize;
            tli = (y-1)*tlsize + 1 : y * tlsize;
            
            regstats(tli,tlj,:,:) = tileco2stat(:,:,:,:,tl);
            kmstat(tli,tlj,:,:) = tilekmstat(:,:,:,:,tl);
            kanmestat(tli,tlj,:) = tilekeranmestat(:,:,:,tl);
            kalmostat(tli,tlj,:) = tilekeralmostat(:,:,:,tl);
            kerdev(tli,tlj,:,:,:) = tilekerdev(:,:,:,:,:,tl);
            kerinstd(tli,tlj,:,:) = tilekerinstd(:,:,:,:,tl);
                        
            clear x y tli tlj
        end
        
        statsubdatafname = strcat(regoutputfiles,"StatsSubData_",num2str(rr),".mat");
        save(statsubdatafname,'regstats','kmstat','kanmestat','kalmostat',...
            'kerdev','kerinstd','desiredstatistics','kernels','pathwayslandcover');
                    
        strcat("done with statistics calculation for subregion ",num2str(rr),...
            " (",num2str(sum(blocklandmask,'all'))," points)")
        

        clear(regvars{:},kervars{:},'blocklandmask','regi','regj')
        clear smallkernels regco2 regstats kmstat kanmestat kalmostat kerdev kerinstd...
            tileco2e tilekernels tileco2stat tilekmstat tilekeranmestat tilekeralmostat ...
            tilekerdev tilekerinstd tilelandmask x y tli tlj tl regi regj pp
               
    end     % end if region has land data
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C. Print results
% *****************

% 1. Recombine maps
% -----------------
mapnames = strings(nlcc,1);
for ll = 1 : nlcc
    map = statisticsmaps;
    
    for rr = 1 : nblocks
        if presland(rr) == true
            osfname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
            load(osfname,'regi','regj');
            subdatafname = strcat(regoutputfiles,"StatsSubData_",num2str(rr),".mat");
            load(subdatafname,'regstats');
            map(regi,regj,:) = regstats(:,:,:,ll);            
        end
    end
    mapnames(ll) = strcat(lower(pathwayslandcover(ll,1)),"2",lower(pathwayslandcover(ll,2)));
    eval(strcat(mapnames(ll)," = map;"))
end
finstafname = strcat(resultslocalfolder,"StatsData.mat");
save(finstafname,mapnames{:},'desiredstatistics','kernels','-v7.3')


% 2. Export CO2e data in geotiff files
% ------------------------------------
lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[lons,lats] = meshgrid(lon05,lat05);
latA = -60;
Antarctica = lats <= latA;

statnames = replace(desiredstatistics," ","_");

R = georasterref('RasterSize',[nlat,nlon],'RasterInterpretation','cells','ColumnsStartFrom',...
    'north','LatitudeLimits',[-90,90],'LongitudeLimits',[-180,180]);
cmptag.Compression = 'LZW';             % LZW compression is better than the default ...    

for ll = 1 : nlcc
    ppi = pathwayslandcover(ll,3);
    pname = pwnames(str2double(ppi));
    for ss = 1 : nstats
        statname = statnames(ss);
        eval(strcat("data = ",mapnames(ll),"(:,:,ss);"))
        data(Antarctica) = NaN; %#ok<SAGROW>
        data(landmask==0) = NaN; %#ok<SAGROW> % just in case
        fname = strcat(resultslocalfolder,pname,"\",upper(mapnames(ll)),"_",statname,".tif");
        geotiffwrite(fname,data,R,'TiffTags',cmptag);
    end
end
clear pname ii ll z data fname
clear(mapnames{:})


% 3. Print figures for statistics on kernels
% ------------------------------------------
% a. build variables
annualstat = NaN(nlat,nlon,nstats);
monthlystat = NaN(nlat,nlon,nmonths,nstats);
for rr = 1 : nblocks
    if presland(rr) == true
        osfname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(osfname,'regi','regj');
        subdatafname = strcat(regoutputfiles,"StatsSubData_",num2str(rr),".mat");
        load(subdatafname,'kmstat','kanmestat');
        annualstat(regi,regj,:) = kanmestat;
        monthlystat(regi,regj,:,:) = kmstat;
    end
end

% b. Get coastines (file and variable names changed with Matlab R2021)
vv = str2double(erase(extract(version,"R" + digitsPattern(4)),"R"));
if vv > 2020
    load coastlines
else
    load coast
    coastlat = lat;
    coastlon = long;
end
hotmap = flipud(colormap(hot));

% b. Statistics on annual means
for ss = 1 : nstats
    statname = statnames(ss);
    
    figure(ss);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    nmap = annualstat(:,:,ss) ./ 100;
    pcolorm(lats,lons,nmap)
    ti = min(nmap,[],'all');
    ta = max(nmap,[],'all');
    if ti < 0
        colormap(hot)
        caxis([round(ti/.5)*.5 0])
    else
        colormap(hotmap)
        caxis([0 round(ta/.5)*.5])
    end
    c = colorbar;
    c.Label.String = {strcat("Annual ",statname," of radiative forcing [W m^-^2]");...
        "corresponding to a 0.01 increase in surface albedo"};
    c.FontSize = 20;
    ax = gca;
    ax.Title.FontSize = 24;
    ax.Title.String = {strcat("Kernels agreements - annual ",statname,...
        " (6 kernels) in Radiative Forcing");strcat("\rm Statistics run on kernels' annual mean ",...
        "values - range is [",num2str(round(ti,2)),",",num2str(round(ta,2)),"]")};
    gridm
    plotm(lat,long,'Color','Black')
    fname = strcat(resultslocalfolder,"Figures\Kernels\KernelAnnualStats_",statname);
    print(strcat(fname,".jpg"),"-djpeg")
    
end


% c. Statistics on monthly values
for ss = 1 : nstats
    statname = statnames(ss);
    
    for m = 1 : nmonths
        figure(ss);clf
        h = gcf;
        h.Units = 'pixels';
        h.Position = [1,31,2560,1333];
        axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
        nmap = monthlystat(:,:,m,ss) ./ 100;
        pcolorm(lats,lons,nmap)
        ti = min(nmap,[],'all');
        ta = max(nmap,[],'all');
        if floor(m/3)*3 == m
            c = colorbar;
            if ti < 0
                colormap(hot)
                caxis([round(ti/.5)*.5 0])
            else
                colormap(hotmap)
                caxis([0 round(ta/.5)*.5])
            end
            c.Label.String = strcat(statname," RF [W m^-^2] / 0.01 \Delta\alpha");
            c.FontSize = 38;
        end
        ax = gca;
        ax.Title.FontSize = 50;
        ax.Title.String = strcat("Kernels - monthly ",statname," - ",monthln(m));
        plotm(coastlat,coastlon,'Color','Black')
        fname = strcat(resultslocalfolder,"Figures\Kernels\KernelMonthlyStats_",...
            statname,".",monthln(m));
        print(strcat(fname,".jpg"),"-djpeg")
    end
    
end


% 4. Table of Kernel "outliers"
% ----------------------------




% Chad Greene (2021). landmask (https://www.mathworks.com/matlabcentral/fileexchange/48661-landmask), MATLAB Central File Exchange. Retrieved October 21, 2021.