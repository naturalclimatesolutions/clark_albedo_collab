% Statistics verifications

% Check a few random points to see if my calculated statistics maps correspond to the same
%   statistics perfomed on the individual kernel maps


clearvars


% 1. User-defined inputs
% ----------------------
computer = "GEO-005199";                % where the routine is run (for directory paths):
desiredstatistics = ...                 % list of statistics to be performed
    ["min","max","mean","median","range","standard deviation"];
nstat = length(desiredstatistics);
npts = 1000;                            % number of random land points on which to verif the stats


% 2. Define inputs
% ----------------
switch computer
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
end

glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'IGBPBiomes','biomes_igbp','latlonscale','nlat','nlon','nlcc','nkernels',...
    'pathwayslandcover','pwnames','kernels','resultslocalfolder','landmask','pwid')

lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[lons,lats] = meshgrid(lon05,lat05);
latA = -60;
Antarctica = lats <= latA;

landmask(Antarctica) = false;
[lat,lon] = find(landmask);


% 3. Calculate/compare statistics from geotiff and statistic maps (on some random land points)
% ---------------------------------------------------------------
for ll = 1 : nlcc
    name = strcat(pathwayslandcover(ll,1),"2",pathwayslandcover(ll,2));
    lc1 = biomes_igbp == pathwayslandcover(ll,1);
    lc2 = biomes_igbp == pathwayslandcover(ll,2);

    folder = pwnames(pwid(ll));
    
    pt = random('discrete uniform',length(lat),[npts,1]);

    kerndata = NaN(npts,nkernels);
    statdata = NaN(npts,nstat);
    
    for kk = 1 : nkernels
        val = readgeoraster(strcat(resultslocalfolder,folder,"\",name,"_",kernels(kk),".tif"));
        for pp = 1 : npts
            kerndata(pp,kk) = val(lat(pt(pp)),lon(pt(pp)));
        end
    end

    for ss = 1 : nstat
        statname = replace(desiredstatistics(ss)," ","_");
        val = readgeoraster(strcat(resultslocalfolder,folder,"\",name,"_",statname,'.tif'));
        for pp = 1 : npts
            statdata(pp,ss) = val(lat(pt(pp)),lon(pt(pp)));
        end
    end

    for ss = 1 : nstat
        stat = desiredstatistics(ss);
        msgerr = strcat("My statistics' calculations don't match for ",stat," in ",name);
        switch stat
            case "min"
                val = min(kerndata,[],2,'omitnan');
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
            case "max"
                val = max(kerndata,[],2,'omitnan');
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
            case "median"
                val = median(kerndata,2,'omitnan');
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
            case "mean"
                val = mean(kerndata,2,'omitnan');
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
            case "standard deviation"
                val = std(kerndata,0,2,'omitnan');
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
            case "range"
                val = range(kerndata,2);
                if isempty(find(val-statdata(:,ss),1)) == 0, error(msgerr); end
        end
    end

end

