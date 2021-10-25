% Stat Verif

% Check a few random points to see if my statistics are correct


clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

computer = "GEO-007264";                % where the routine is run (for directory paths):
                                        %   GEO-007264 (my laptop)
                                        %   GEO-005199 (old supercomputer)
                                        %   GEO-006381 (new supercomputer)
desiredstatistics = ...                 % list of statistics to be performed
    ["min","max","mean","median","range","standard deviation"];
nstat = length(statistics);
npts = 1000;

% **********************


% 2. Define inputs
% ----------------
switch computer
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
        mapdir = strcat('G:\GlobalAlbedo\');
        resultdir = rootdir;
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
        mapdir = strcat('G:\TNC\GlobalAlbedo\');
        resultdir = 'E:\NCS-GlobalAlbedo\';
end

glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'IGBPBiomes','biomes_igbp','latlonscale','nblocks','nlat','nlon','nlcc',...
    'nmonths','nkernels','npw','presland','pathwayslandcover','pwnames','pwid','kernels',...
    'regoutputfiles','resultslocalfolder','landmask','blocksize','kernelscale')


[lat,lon] = find(landmask);

pt = random('uniform',1,length(lat),[npts,1]);

kerndata = NaN(npts,nkernels);
for kk = 1 : nkernels
    fname = strcat(resultslocalfolder,'Forest2Grass\EBF2GRA_',kernels(kk),'.tif');
    val = readgeoraster(fname);
    kerndata(:,kk) = var(lat(pt),lon(pt));
end

statdata = NaN(npts,nstat);
for ss = 1 : nstat
    statname = replace(desiredstatistics(ss)," ","_");
    fname = strcat(resultslocalfolder,'Forest2Grass\EBF2GRA_',statname(ss),'.tif');
    val = readgeoraster(fname);
    statdata(:,ss) = var(lat(pt),lon(pt));
end


for ss = 1 : nstat
    stat = desiredstatistics(ss);
    switch stat
        case "min"
            val = min(kerndata,[],2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
        case "max"
            val = max(kerndata,[],2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
        case "median"
            val = median(kerndata,2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
        case "mean"
            val = mean(kerndata,2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
        case "standard deviation"
            val = std(kerndata,0,2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
        case "range"
            val = range(kerndata,2,'omitnan');
            if isempty(val-statdata(:,ss)) == 0
                warning(strcat("My statistics calculation don't match for ",stat));
            end
    end
end



    


