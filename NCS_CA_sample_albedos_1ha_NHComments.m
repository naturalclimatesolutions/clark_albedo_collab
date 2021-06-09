% Step One: Generate Climate Data Structure w/ Albedos, Snow, Radiation, and Kernels

clarvars

% NH : DEFINE ARRAY SIZE
latlonscale = 0.05;
lats_pt05deg = (-90.0+0.025):0.05:(90-0.025);
lons_pt05deg = (-180.0+0.025):0.05:(180-0.025);
nlat = length(lats_pt05deg);
nlon = length(lons_pt05deg);
lats2dhi = repmat(lats_pt05deg',1,nlon);
lons2dhi = repmat(lons_pt05deg,nlat,1);

% lat/lon Canada
latmax = 84;   latmin = 41;
lonmin = -142; lonmax = -52;
latscanada = find(lats_pt05deg < latmin, 1, 'last' ) : find(lats_pt05deg > latmax, 1 , 'first' );
lonscanada = find(lons_pt05deg < lonmin, 1, 'last' ) : find(lons_pt05deg > lonmax, 1 , 'first' );


% ALBEDO FILE INFO
fillvalue_albedo = 32767;
scale_albedo = 0.001;
biomes_albedo = {'WAT','ENF','EBF','DNF','DBF','MF','CSH','OSH','WSA','SAV','GRA','WET','CRO','URB','MOS','SNO','BAR'};
nbiomes = length(biomes_albedo);
monthnames = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
lcmonthnames = lower(monthnames);
nmonths = length(monthnames);       %% NH - I like have variabels instead of hardcoved numbers
% NH - added this
albedo_types = ["snow-free black-sky","snow-free white-sky",...
    "snow-covered black-sky","snow-covered white-sky"];
nalb = length(albedo_types);
% paths to extracted data files
filepathdatasnwcvrd = 'D:\DATA\ALBEDO\ALBEDO_ATLAS\EXTRACTED\snow_covered_hierarchical_v3\';
filepathdatasnwfree = 'D:\DATA\ALBEDO\ALBEDO_ATLAS\EXTRACTED\snow_free_hierarchical_v3\';
% albedo header files should always be the same - so we open just one
filepathheadersnwcvrd = 'D:\DATA\ALBEDO\ALBEDO_ATLAS\EXTRACTED\snow_covered_hierarchical_v3\';
fnameheader = 'hierarchical_snow_albedo_igbp_0.05.bsa_shortwave.Jan.hdr'; 
headerfile_albedo = [filepathheadersnwcvrd,fnameheader];


% NH - DEFINE STORAGE ARRAYS
snowcover_all = zeros(nlat,nlon,nmonths);
direct_solar = zeros(180,360,nmonths);
diffuse_solar = zeros(180,360,nmonths);
solar_direct_hires  = zeros(nlat,nlon,nmonths);
solar_diffuse_hires = zeros(nlat,nlon,nmonths);
frac_diff_hires     = zeros(nlat,nlon,nmonths);
albkern_hires = zeros(nlat,nlon,nmonths);
albedosout = NaN(nlat,nlon,nmonths,nbiomes,nalb);


% FILTER WITH A LAND MASK FROM MODIS LAND COVER
disp('Reading Land Mask')
path = 'D:\DATA\ALBEDO\LANDCOVER\';
this = 'MCD12C1.A2016001.006.2018324172410.hdf';
fname = [path,this];
S = hdfinfo(fname,'eos');
varname = S.Grid.DataFields.Name;
A = hdfread(fname,varname);
A = flipud(A);
landmask = A;
landmask(A>0) = 1;
clear A

% SNOW COVER DATA
% ADDRESSES SNOW DATASET ISSUES
% 1) fix no data issues above arctic circle, set to 100%    %% NH - October to March (?)
% 2) fix NaN land points, set to 0% snow cover
disp('Reading Snow Cover')
srcdir = 'C:\Users\cwilliams\Documents\workspace\CLARK\TNC\CA_NCS\SNOW_DATA\';
snow_prefix = 'snowcover_climatology_';
% snowmonthstr = {'jan','feb','mar','apr','may',...
%     'jun','jul','aug','sep','oct','nov','dec'};
for m = 1 : nmonths
    fname = [srcdir,snow_prefix,lcmonthnames{m},'.mat'];
    load(fname,'snowcover');    %% NH I would add "snowcover in the read to avoid warnings
    temp = sum(snowcover,2,'omitnan');  %% NH I suspect "snowcover is a [lat,lon] array
    if m < 4 || m > 9
        x = find(temp>0,1,'last') -10; % northernmost non-zero value minus 5 (NH 10?) to address taper issue
        snowcover(x:nlat,:) = 100;  %% NH only sets snow to 100% from October to March
    end
    snowcover(isnan(snowcover)) = 0;
    snowcover(landmask==0) = NaN;
    snowcover_all(:,:,m) = snowcover./100; % fraction of month-pixel snow covered
end



% SOLAR RADIATION DATA
disp('Reading Solar Radiation Data')
% solarmonthstr= {'jan','feb','mar','apr','may','jun','jul',...
%     'aug','sep','oct','nov','dec'};
for m = 1 : nmonths
    fname = ['C:\Users\cwilliams\Documents\workspace\CLARK\ALBEDO\solarflux\Visible_plus_NIR_beam_downward_',lcmonthnames{m},'.rst'];
    fid = fopen(fname);
    clear A
    A = fread(fid,'single');
    A = reshape(A,360,180);
    A = A'; A = flipud(A);
    fclose(fid);
    direct_solar(:,:,m) = A;
    % figure(1);clf; pcolor(A);shading flat; colorbar
    % titletxt = solarmonthstr(m); title(titletxt); pause
    fname = ['C:\Users\cwilliams\Documents\workspace\CLARK\ALBEDO\solarflux\Visible_plus_NIR_diffuse_downward_',lcmonthnames{m},'.rst'];
    fid = fopen(fname);
    clear A
    A = fread(fid,'single');
    A = reshape(A,360,180);
    A = A'; A = flipud(A);
    fclose(fid);
    diffuse_solar(:,:,m) = A;
    % figure(1);clf; pcolor(A);shading flat; colorbar
    % titletxt = solarmonthstr(m); title(titletxt)
end
frdiff = diffuse_solar ./ (diffuse_solar + direct_solar);


% RESAMPLE RADIATION DATA TO 0.05 DEG GRID TO MATCH ALBEDO AND SNOW DATA
% lats_1deg = (-90.0+0.5):1:(90-0.5);
% lons_1deg = (-180.0+0.5):1:(180-0.5);
% lats2d1deg = repmat(lats_1deg',1,360);
% lons2d1deg = repmat(lons_1deg,180,1);
disp('Resampling Solar Radiation Data')
% lats_pt05deg = (-90.0+0.025):0.05:(90-0.025);     %% NH - I defined these at the start
% lons_pt05deg = (-180.0+0.025):0.05:(180-0.025);
% lats2dhi = repmat(lats_pt05deg',1,7200);  
% lons2dhi = repmat(lons_pt05deg,3600,1);
% solar_direct_hires  = zeros(3600,7200,12);
% solar_diffuse_hires = zeros(3600,7200,12);
% frac_diff_hires     = zeros(3600,7200,12);
scaleratio = 1/0.05;
nx = nlat / scaleratio;
ny = nlon / scaleratio;
for m = 1 : nmonths
    for x = 1 : nx
        for y = 1 : ny
            xdestinations = (x-1) *scaleratio +1 : x *scaleratio;
            ydestinations = (y-1) *scaleratio +1 : y *scaleratio;
            solar_direct_hires(xdestinations,ydestinations,m)  = direct_solar(x,y,m);
            solar_diffuse_hires(xdestinations,ydestinations,m) = diffuse_solar(x,y,m);
            frac_diff_hires(xdestinations,ydestinations,m)     = frdiff(x,y,m);
        end
    end
end



% ALBEDO KERNELS DATA
% Shell et al. 2008 dataset
% W m-2 per 0.01 change in fractional surface albedo
disp('Reading Albedo Kernels')
% path = 'D:\DATA\ALBEDO\RADIATIVE_KERNELS\CAM3\';
% thisfile = 'CAM3_albedo_sw_kernel.nc';
% fname = [path,thisfile];
% % finfo = ncinfo(fname);
% % varinfo = finfo.Variables;
% % disp(varinfo(3))
% albkernorig = ncread(fname,'monkernel');
% albkernorig = permute(albkernorig,[2 1 3]);
% lats_albkern1 = ncread(fname,'lat');
% lons_albkern1 = ncread(fname,'lon');
% stepindegrees = round(mean(diff(lons_albkern1)).*10000)./10000;
% lons_albkern1 = lons_albkern1 - 180 + stepindegrees/2;
% lats2dkern = repmat(lats_albkern1,1,128);
% lons2dkern = repmat(lons_albkern1',64,1);
% %for m =1:12; figure(1);clf;pcolor(lons2dkern,lats2dkern,albkernorig(:,:,m));shading flat;colorbar;title(num2str(m));pause;end


path = 'D:\DATA\ALBEDO\RADIATIVE_KERNELS\CAM5\';
thisfile = 'alb.kernel.nc';
%FSNS : Net Shortwave Flux at Surface, all-sky (positive: down)
%FSNSC: Net Shortwave Flux at Surface, Clear-sky (positive: down)
%FSNT : Net Shortwave Flux at Top-of-atmosphere, all-sky (positive: down)
%FSNTC: Net Shortwave Flux at Top-of-atmosphere, Clear-sky (positive: down)
fname = [path,thisfile];
%finfo = ncinfo(fname);
%varinfo = finfo.Variables;
%disp(varinfo(3))
albkern2 = ncread(fname,'FSNT');
albkernorig2 = permute(albkern2,[2 1 3]);
%wrap data by 180 degrees of longitude
x=size(albkernorig2);
albkernorig2wrapped = cat(2,albkernorig2(:,x(2)/2:x(2),:),albkernorig2(:,1:x(2)/2-1,:));
lats_albkern2 = ncread(fname,'lat');
lons_albkern2 = ncread(fname,'lon');
stepindegreesy = round(mean(diff(lons_albkern2)).*10000)./10000;
stepindegreesx = round(mean(diff(lats_albkern2)).*10000)./10000;
lons_albkern2 = lons_albkern2 - 180 + mean(diff(lons_albkern2))/2;
lats2dkern2 = repmat(lats_albkern2,1,288);
lons2dkern2 = repmat(lons_albkern2',192,1);

% RESAMPLE KERNEL DATA TO 0.05 DEG GRID TO MATCH ALBEDO AND SNOW DATA
disp('Resampling Albedo Kernels')
% albkern1_hires = zeros(3600,7200,12);
scaleratioy = stepindegreesy / 0.05;
scaleratiox = stepindegreesx / 0.05;
nx = nlat / scaleratiox;
ny = nlon / scaleratioy;
for m = 1 : nmonths
    for x = 1 : nx
        for y = 1 : ny
            xdestinations = floor((x-1).*scaleratiox+1):ceil(x.*scaleratiox);
            ydestinations = floor((y-1).*scaleratioy+1):ceil(y.*scaleratioy);
            albkern_hires(xdestinations,ydestinations,m)  = albkernorig2wrapped(x,y,m);
        end
    end
end

% clip snow, radiation, and kernels to Canada
disp('Clipping Files to Canada')
% latmax = 84;   latmin = 41;
% lonmin = -142; lonmax = -52;
% latscanada = max(find(lats_pt05deg < latmin)):min(find(lats_pt05deg > latmax));
% lonscanada = max(find(lons_pt05deg < lonmin)):min(find(lons_pt05deg > lonmax));
solardir = solar_direct_hires(latscanada,lonscanada,:);
solardif = solar_diffuse_hires(latscanada,lonscanada,:);
fracdiff = frac_diff_hires(latscanada,lonscanada,:);
snowcover = snowcover_all(latscanada,lonscanada,:);
albkern = albkern_hires(latscanada,lonscanada,:);
landmask = landmask(latscanada,lonscanada,:);
lats2dcanada = lats2dhi(latscanada,lonscanada);
lons2dcanada = lons2dhi(latscanada,lonscanada);
clear solar_direct_hires solar_diffuse_hires frac_diff_hires snowcover_all lats2dhi lons2dhi 
clear frdiff diffuse_solar direct_solar

% for m=1:12;pcolor(albkern(:,:,m));shading flat;colorbar;pause;end

% Not needed now that we are performing this on a 1 ha grid
% % load ecoregions file
% % File Contains: allGrid allRefVec lats2dcanada lons2dcanada
% load('D:\DATA\ALBEDO\ECOREGIONS\CanadaEcoregions','allGrid','allRefVec')
% ecoregionnums = unique(allGrid); % should be 194
% necoregions = size(ecoregionnums,1);

% loop over months to sample the datasets
for m = 1:12 % months of year
    
    disp(['Month = ',num2str(m)])
    
    % It looks like we will be able to hold the snow and radiation datasets in
    % memory while sampling the albedo data
    thissnowcover = snowcover(:,:,m);
    thissolardir  = solardir(:,:,m);
    thissolardif  = solardif(:,:,m);
    thisalbkern   = albkern(:,:,m);
    
    % read in this month's albedo datasets    
    mmm = monthnames{m};
    
    % snow-free black-sky
    disp('Reading Albedos')
    fname = ['hierarchical_albedo_igbp_0.05.bsa_shortwave.',mmm,'.mean'];
    datafile = [filepathdatasnwfree,fname];
    %datafile = 'D:\DATA\ALBEDO\ALBEDO_ATLAS\EXTRACTED\snow_free_hierarchical_v3\hierarchical_albedo_igbp_0.05.bsa_shortwave.Jan.mean';
    [varargout] = enviread(datafile,headerfile_albedo);
    varargout(varargout == fillvalue_albedo) = NaN;
    varargout = varargout.*scale_albedo;
    bsa_nosnow_allbiomes = varargout(latscanada,lonscanada,:);
    clear varargout
    
    % snow-free white-sky
    fname = ['hierarchical_albedo_igbp_0.05.wsa_shortwave.',mmm,'.mean'];
    datafile = [filepathdatasnwfree,fname];
    [varargout] = enviread(datafile,headerfile_albedo);
    varargout(varargout == fillvalue_albedo) = NaN;
    varargout = varargout.*scale_albedo;
    wsa_nosnow_allbiomes = varargout(latscanada,lonscanada,:);
    clear varargout
    
    % snow-covered black-sky
    fname = ['hierarchical_snow_albedo_igbp_0.05.bsa_shortwave.',mmm,'.mean'];
    datafile = [filepathdatasnwcvrd,fname];
    [varargout] = enviread(datafile,headerfile_albedo);
    varargout(varargout == fillvalue_albedo) = NaN;
    varargout = varargout.*scale_albedo;
    bsa_snow_allbiomes = varargout(latscanada,lonscanada,:);
    clear varargout

    % snow-covered white-sky
    fname = ['hierarchical_snow_albedo_igbp_0.05.wsa_shortwave.',mmm,'.mean'];
    datafile = [filepathdatasnwcvrd,fname];
    [varargout] = enviread(datafile,headerfile_albedo);
    varargout(varargout == fillvalue_albedo) = NaN;
    varargout = varargout.*scale_albedo;
    wsa_snow_allbiomes = varargout(latscanada,lonscanada,:);
    clear varargout
    
    disp('Populating Data Structures')
    for b = 1:17
        albedosout(:,:,m,b,1) = bsa_nosnow_allbiomes(:,:,b);
        albedosout(:,:,m,b,2) = wsa_nosnow_allbiomes(:,:,b);
        albedosout(:,:,m,b,3) = bsa_snow_allbiomes(:,:,b);
        albedosout(:,:,m,b,4) = wsa_snow_allbiomes(:,:,b);
    end    
    
    snowout(:,:,m)    = thissnowcover;  %% NH I dont see the point of these ... no diff with snowcover ...
    solarout(:,:,m,1) = thissolardir;
    solarout(:,:,m,2) = thissolardif;
    kernelsout(:,:,m) = thisalbkern;

end % on month

stophere

% Gap-filling the albedo dataset
% if one of the albedo types is missing values, assign an alternate
% from adjacent months of a similar condition
for b = 1:17
    for a = 1:4
        for m = 1:12
            
            % first time through look to adjacent months
            albedosfill = squeeze(albedosout(:,:,m,b,a));
            bad = find(isnan(albedosfill) == 1);
            mfillindex = [m-1,m+1];
            if m == 1; mfillindex = [2,12];end
            if m ==12; mfillindex = [1,11];end
            fillvals = squeeze(nanmean(albedosout(:,:,mfillindex,b,a),3));
            albedosfill(bad) = fillvals(bad);
            albedosout(:,:,m,b,a) = albedosfill;
            
            % second time through look to the 2 adjacent months
            albedosfill = squeeze(albedosout(:,:,m,b,a));
            bad = find(isnan(albedosfill) == 1);
            mfillindex = [m-2,m-1,m+1,m+2];
            if m == 1; mfillindex = [2,3,11,12];end
            if m == 2; mfillindex = [1,3,4,12];end
            if m ==11; mfillindex = [1,9,10,12];end
            if m ==12; mfillindex = [1,2,10,11];end
            fillvals = squeeze(nanmean(albedosout(:,:,mfillindex,b,a),3));
            albedosfill(bad) = fillvals(bad);
            albedosout(:,:,m,b,a) = albedosfill;
            
            % third time through look to the 3 adjacent months
            albedosfill = squeeze(albedosout(:,:,m,b,a));
            bad = find(isnan(albedosfill) == 1);
            mfillindex = [m-3,m-2,m-1,m+1,m+2,m+3];
            if m == 1; mfillindex = [2,3,4,10,11,12];end
            if m == 2; mfillindex = [1,3,4,5,11,12];end
            if m == 3; mfillindex = [1,2,4,5,6,12];end
            if m ==10; mfillindex = [1,7,8,9,11,12];end
            if m ==11; mfillindex = [1,2,8,9,10,12];end
            if m ==12; mfillindex = [1,2,3,9,10,11];end
            fillvals = squeeze(nanmean(albedosout(:,:,mfillindex,b,a),3));
            albedosfill(bad) = fillvals(bad);
            albedosout(:,:,m,b,a) = albedosfill;
        end
    end
end

fnameout = 'Albedo_Climate_Data_Canada_Gridded.mat';
save(fnameout,'albedosout','snowout','solarout','kernelsout','-v7.3');

% 
% % Check that the dataarrays match the source data
% for m=1:3%:12
%     checkdata = snowout;
%     sourcedata = snowcover(:,:,m);
%     sourcedata(1,1) = 0;
%     datacheck = snowout(:,:,m);
%     datacheck(1,1) = 0;
%     figure(10);clf
%     subplot(2,1,1)
%     pcolor(sourcedata);shading flat;colorbar
%     subplot(2,1,2)
%     pcolor(datacheck);shading flat;colorbar
%     pause
% end

% thisbiome = strmatch('DNF',biomes_albedo);
% for m=1:3%:12
%     mmm = monthnames{m};
%     fname = ['hierarchical_snow_albedo_igbp_0.05.bsa_shortwave.',mmm,'.mean'];
%     datafile = [filepathdatasnwcvrd,fname];
%     [varargout] = enviread(datafile,headerfile_albedo);
%     varargout(varargout == fillvalue_albedo) = NaN;
%     varargout = varargout.*scale_albedo;
%     bsa_snow_allbiomes = varargout(latscanada,lonscanada,:);
%     sourcedata = bsa_snow_allbiomes(:,:,thisbiome);
%     clear varargout
%     %checkdata = squeeze(albedosout(:,:,thisbiome,3)); % 4th dim = 1 is bsa nosnow; 3 is bsa snow covered
%     sourcedata(1,1) = 0;
%     datacheck = squeeze(albedosout(:,:,m,thisbiome,3));
%     datacheck(1,1) = 0;
%     figure(10);clf
%     subplot(2,1,1)
%     pcolor(sourcedata);shading flat;colorbar
%     subplot(2,1,2)
%     pcolor(datacheck);shading flat;colorbar
%     pause
% end
