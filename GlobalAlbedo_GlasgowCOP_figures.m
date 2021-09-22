% NCS - Global Albe/do - TNC-Clark Collaboration - Glasgow COP figures
%
% This code is used to produce figures of all inputs and outputs used/produced for the Glasgow COP
%   meeting of November 2021.
%
% Created 9/15/2021 by Natalia Hasler, Clark University

% Last modified: NH 9/22/2021
%                (See Modification history at end of file)



clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Preliminary (get general data)
% ***************

computer = "GEO-005199";                % change this one if needed

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
load(glodatafname,'albsn','IGBPBiomes','biomes_igbp','latlonscale','nblocks','nlat','nlon','nlcc',...
    'nmonths','nkernels','npw','presland','pathwayslandcover','pwnames','pwid','kernels',...
    'regoutputfiles','resultslocalfolder','landmask')

monthln = ["January","February","March","April","May","June","July","August",...
    "September","October","November","December"];

lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[lons,lats] = meshgrid(lon05,lat05);

latA = -60;
Antarctica = lats <= latA;

greyval = [0.7,0.7,0.7];
lightgreyval = [.9,.9,.9];
jetmap = colormap(jet);
lastblue = find(jetmap(:,1)==0,1,'last');
firstyellow = find(jetmap(:,1)==1,1,'first');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B. Input data Analysis
% ***********************

% 1. Plot SNOW COVER
% ------------------
% a. Get data
monthlymap =  NaN(nlat,nlon,nmonths);
map = NaN(nlat,nlon);
load(glodatafname,'MisSnow');

snowbar = flipud(jetmap(1:lastblue,:));
snowbar = cat(1,greyval,snowbar);

for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname,'snowcover','regi','regj');
        
        for m = 1 : nmonths
            monthlymap(regi,regj,m) = snowcover(:,:,m);
        end
        map(regi,regj) = mean(snowcover,3,'omitnan');
        
        clear snowcover regi regj m
    end
end

% b. Print annual map
figure(1);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,31,2560,1333];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
load coast
pcolorm(lats,lons,map)
colormap(snowbar)
c = colorbar;
c.Label.String = 'Annual mean snow cover fraction';
c.FontSize = 20;
ax = gca;
ax.Title.FontSize = 24;
ax.Title.String = "Annual mean snow cover fraction";
gridm
plotm(lat,long,'Color','Black')
fname = strcat(resultslocalfolder,"Figures\AnnualSnow");
print(strcat(fname,".jpg"),"-djpeg")

% c. Print Monthly maps
for m = 1 : nmonths
    figure(10+m);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','off','ParallelLabel','off');
    pcolorm(lats,lons,monthlymap(:,:,m))
    caxis([0 1])
    colormap(snowbar)
    if floor(m/3)*3 == m
        c = colorbar;
        c.Label.String = 'Monthly mean snow cover fraction';
        c.FontSize = 38;
    end
    plotm(lat,long,'Color','Black')
    ax = gca;
    ax.Title.FontSize = 50;
    ax.Title.String = monthln(m);
    fname = strcat(resultslocalfolder,"Figures\MonthlySnow.",monthln(m));
    print(strcat(fname,".jpg"),"-djpeg")
end

% d. Plot map of missing data
misdata = sum(double(MisSnow),3);
misdata(landmask==0) = NaN;
x = flipud(colormap(hot));
mint = floor(length(x)/11);
misbar = x(1:mint:length(x),:);
if length(misbar) > 11; misbar = misbar(1:length(misbar)-1,:); end
misbar = cat(1,lightgreyval,misbar);

figure(2);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,31,2560,1333];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
pcolorm(lats,lons,misdata)
colormap(misbar)
c = colorbar;
c.Label.String = 'Number of months of missing data';
c.FontSize = 20;
ax = gca;
ax.Title.FontSize = 24;
ax.Title.String = "Snow cover missing data";
gridm
fname = strcat(resultslocalfolder,"Figures\SnowMissingData");
print(strcat(fname,".jpg"),"-djpeg")

clear monthlymap map snowbar x misbar fname c ax h mint misdata lastblue rr
close all



% 2. Plot WEIGHTED ALBEDO
% -----------------------
% a. Prepare map scale and colorbar
albedovars = cellstr(albsn);
albedoscale = [0:0.05:0.4,0.5,0.7,1];
nbcc = length(albedoscale)-1;
x = colormap(hot);
mapint = floor(length(x) / nbcc);
cathot = x(1:mapint:length(x),:);
if length(cathot) > nbcc, cathot = cathot([1:length(cathot)-2,length(cathot)],:); end

% b. Loop through biome-type
for ff = 1 : nbiomes
    monthlymap =  NaN(nlat,nlon,nmonths);
    map = NaN(nlat,nlon);
    
    for rr = 1 : nblocks
        if presland(rr) == true
            subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
            load(subdatafname,albedovars{:},'snowcover','diffusefraction',...
                'regi','regj','blocklandmask');
            
            % c. calculate monthly and mean annual weighted albedo values
            for m = 1 : nmonths
                snow = snowcover(:,:,m);
                temp = snow;
                temp(blocklandmask==0) = 100;
                if isempty(find(isnan(temp),1)) == 0
                    error(strcat("I still have missing values of snow in ",monthln(m),...
                        " and subregion ",num2str(rr)));
                end
                clear temp
                Wda = squeeze(...
                    nosnowblacksky(:,:,m,ff) .* (1 - snow) .* (1 - diffusefraction(:,:,m)) + ...
                    nosnowwhitesky(:,:,m,ff) .* (1 - snow) .* diffusefraction(:,:,m) + ...
                    snowcovblacksky(:,:,m,ff) .* snow .* (1 - diffusefraction(:,:,m)) + ...
                    snowcovwhitesky(:,:,m,ff) .* snow .* diffusefraction(:,:,m));
                monthlymap(regi,regj,m) = Wda;
            end
            map(regi,regj) = mean(Wda,3,'omitnan');
            
            clear(albedovars{:})
            clear snowcover diffusefraction regi regj m snow
        end
    end
    
    % d. Print annual map
    figure(ff);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    catmap = discretize(map,albedoscale);
    val1 = round(min(map,[],'all'),2);
    val2 = round(max(map,[],'all'),2);    
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    load coast
    pcolorm(lats,lons,catmap)
    caxis([1 nbcc+1])
    colormap(cathot)
    c = colorbar;
    c.Ticks = 1:nbcc+1;
    c.TickLabels = num2cell(albedoscale);
    c.Label.String = 'Annual mean weighted blue-sky albedo';
    c.FontSize = 20;
    ax = gca;
    ax.Title.FontSize = 24;
    ax.Title.String = {strcat(IGBPBiomes(ff)," - mean annual weighted blue-sky albedo values");...
        strcat("\rm range is [",num2str(val1),",",num2str(val2),"]")};
    gridm
    plotm(lat,long,'Color','Black')
    fname = strcat(resultslocalfolder,"Figures\AnnualAlbedo_",biomes_igbp(ff));
    print(strcat(fname,".jpg"),"-djpeg")

    %. e. Print monthly map
    for m = 1 : nmonths
        figure(100+m);clf
        h = gcf;
        h.Units = 'pixels';
        h.Position = [1,31,2560,1333];
        catmap = discretize(monthlymap(:,:,m),albedoscale);
        val1 = round(min(monthlymap(:,:,m),[],'all'),2);
        val2 = round(max(monthlymap(:,:,m),[],'all'),2);
        axesm('MapProjection','robinson','Frame','on','MeridianLabel','off','ParallelLabel','off');
        pcolorm(lats,lons,catmap)
        caxis([1 nbcc+1])
        colormap(cathot)
        if floor(m/3)*3 == m
            c = colorbar;
            c.Ticks = 1:nbcc+1;
            c.TickLabels = num2cell(albedoscale);
            c.Label.String = 'Monthly weighted blue-sky albedo';
            c.FontSize = 38;
        end
        if biomes_igbp(ff) == "SNO" || biomes_igbp(ff) == "WAT", plotm(lat,long,'Color','Black'); end
        ax = gca;
        ax.Title.FontSize = 50;
        ax.Title.String = strcat(monthln(m)," [",num2str(val1),",",num2str(val2),"]");
        fname = strcat(resultslocalfolder,"Figures\MonthlyAlbedo_",biomes_igbp(ff),".",monthln(m));
        print(strcat(fname,".jpg"),"-djpeg")
    end
    
    pause
    close all
  
end

clear monthlymap albedovars albedoscale nbcc x mapint cathot fname c ax h rr m catmap val1 val2



% 3. Plot KERNELS
% ---------------
% a. Print mean original data
%    (this is not exactly right, because given that lat/lon seem to be centers, I would need to
%    add/modify the lat/lon to match expected corner values for Matlab projections)
load(glodatafname,'kernelfilenames','kernelvarnames')
radiativekernel = strcat(mapdir,"RadiativeKernels\");

for kk = 1 : nkernels
    subdatafname = strcat(radiativekernel,kernels(kk),"\",kernelfilenames(kk));
    korig = permute(ncread(subdatafname,kernelvarnames(kk)),[2 1 3]);
    kolat = double(ncread(subdatafname,'lat'));
    kolon = double(ncread(subdatafname,'lon'));
    [klons,klats] = meshgrid(kolon,kolat);
    mker = mean(korig,3);
    figure(kk);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','off','ParallelLabel','off');
    pcolorm(klats,klons,mker)
    caxis([-3 0])
    colormap(hot)
    c = colorbar;
    c.Label.String = 'Radiative forcing [W m^-^2] for a 0.01 albedo increase';
    c.FontSize = 20;
    ax = gca;
    ax.Title.FontSize = 24;
    ax.Title.String = strcat(kernels(kk)," - Original dataset radiative forcing");
    gridm
    plotm(lat,long,'Color','Black')
    fname = strcat(resultslocalfolder,"Figures\OriginalKernel_",kernels(kk));
    print(strcat(fname,".jpg"),"-djpeg")
    
    clear subdatafname korig kolat kolon klons klats mker h c ax fname
end
clear radiativekernel kernelfilenames kernelvarnames kk

% b. Calculate annual and monthly ratios between kernels
monthlymap =  NaN(nlat,nlon,nmonths,nkernels);
map = NaN(nlat,nlon,nkernels);
kvn = cellstr(lower(kernels));

for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname,kvn{:},'regi','regj');
        
        for kk = 1 : nkernels
            eval(strcat("kname = ",lower(kernels(kk)),";"));
            for m = 1 : nmonths
                monthlymap(regi,regj,m,kk) = kname(:,:,m);
            end
            map(regi,regj,kk) = mean(kname,3,'omitnan');
        end
        
        clear(kvn{:})
        clear kk kname regi regj m
    end
end

% c. Print annual maps
rvals = [0 1/4 1/3 1/2 1/1.5 .95 1.05 1.5 2 3 4 5];
nbc = length(rvals) - 1;
if floor(nbc/2) == nbc/2
    error("I need to have an odd number of bins here to be simetrical with a zero bin");
end
nbcc = (nbc-1)/2;
blueint = floor(lastblue/nbcc);
bluescale = jetmap(1:blueint:lastblue,:);
nlb = numel(bluescale)/3;
if nlb > nbcc, bluescale = bluescale([1:nlb-2,nlb],:); end
yellowint = floor((length(jetmap)-firstyellow)/nbcc);
yellowscale = jetmap(firstyellow:yellowint:length(jetmap),:);
nly = numel(yellowscale)/3;
if nly > nbcc, yellowscale = yellowscale([1,3:nly],:); end
catblind = cat(1,bluescale,greyval,yellowscale);

for kk = 1 : nkernels
    figure(kk);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    nmap = map(:,:,kk);
    nmap(landmask==0) = NaN;
    pcolorm(lats,lons,nmap)
    caxis([-3 0])
    colormap(hot)
    c = colorbar;
    c.Label.String = 'Radiative forcing [W m^-^2] for a 0.01 surface albedo increase';
    c.FontSize = 20;
    ax = gca;
    ax.Title.FontSize = 24;
    ax.Title.String = strcat(kernels(kk)," - Annual TOA (resampled) radiative forcing ",...
        "due to changes in surface albedo");
    gridm
    plotm(lat,long,'Color','Black')
    fname = strcat(resultslocalfolder,"Figures\ResampledKernel_",kernels(kk));
    print(strcat(fname,".jpg"),"-djpeg")
    
    if kk == 2, continue; end
    rmap = nmap ./ map(:,:,2);
    rmap(landmask==0) = NaN;
    val1 = round(min(rmap,[],'all'),2);
    val2 = round(max(rmap,[],'all'),2);
    rmcat = discretize(rmap,rvals);
    
    figure(10+kk);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,31,2560,1333];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    pcolorm(lats,lons,rmcat)
    caxis([1 length(rvals)])
    colormap(catblind)
    c = colorbar;
    c.Ticks = 1:length(rvals);
    c.TickLabels = num2cell(rvals);
    c.FontSize = 20;
    c.Label.String = 'radiative forcing ratio';
    c.Ticks = 1:length(rvals);
    c.TickLabels = num2cell(round(rvals,2));
    plotm(lat,long,'Color','Black')
    ax = gca;
    ax.Title.FontSize = 24;
    ax.Title.String = strcat("Annual mean TOA radiative forcing ratio ",kernels(kk)," / CAM5",...
        " [",num2str(val1),",",num2str(val2),"]");
    gridm
    plotm(lat,long,'Color','Black')
    fname = strcat(resultslocalfolder,"Figures\KernelRatio_",kernels(kk));
    print(strcat(fname,".jpg"),"-djpeg")
  
end    
    
% d. Print monthly maps
for kk = 1 : nkernels
    for m = 1 : nmonths
        figure(m);clf
        h = gcf;
        h.Units = 'pixels';
        h.Position = [1,31,2560,1333];
        axesm('MapProjection','robinson','Frame','on','MeridianLabel','off','ParallelLabel','off');
        nmap = monthlymap(:,:,m,kk);
        nmap(landmask==0) = NaN;
        pcolorm(lats,lons,nmap)
        caxis([-5 0])
        colormap(hot)
        if floor(m/3)*3 == m
            c = colorbar;
            c.Label.String = 'Radiative forcing [W m^-^2] / 0.01 \Delta\alpha';
            c.FontSize = 38;
        end
        ax = gca;
        ax.Title.FontSize = 50;
        ax.Title.String = strcat(kernels(kk)," - ",monthln(m));
        plotm(lat,long,'Color','Black')
        fname = strcat(resultslocalfolder,"Figures\ResampledKernel.",kernels(kk),".",monthln(m));
        print(strcat(fname,".jpg"),"-djpeg")
        
        if kk == 2, continue; end
        nmap2 = nmap;
        nmap(nmap>-10^-1) = -10^-1;
        cam5m = monthlymap(:,:,m,2);
        cam5m(cam5m>-10^-1) = -10^-1;
        rmap = nmap ./ cam5m;
        rmap(landmask==0) = NaN;
        val1 = round(min(rmap,[],'all'),2);
        val2 = round(max(rmap,[],'all'),2);
        rmcat = discretize(rmap,rvals);
        
        figure(20+m);clf
        h = gcf;
        h.Units = 'pixels';
        h.Position = [1,31,2560,1333];
        axesm('MapProjection','robinson','Frame','on','MeridianLabel','off','ParallelLabel','off');
        pcolorm(lats,lons,rmcat)
        caxis([1 length(rvals)])
        colormap(catblind)
        if floor(m/3)*3 == m
            c = colorbar;
            c.Ticks = 1:length(rvals);
            c.TickLabels = num2cell(round(rvals,2));
            c.Label.String = 'Radiative forcing ratio';
            c.FontSize = 38;
        end
        ax = gca;
        ax.Title.FontSize = 50;
        ax.Title.String = strcat(monthln(m)," [",num2str(val1),",",num2str(val2),"]");
        fname = strcat(resultslocalfolder,"Figures\KernelRatio.",kernels(kk),".",monthln(m));
        print(strcat(fname,".jpg"),"-djpeg")
    end  
end    






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C. Results in CO2e
% *******************
datafname = strcat(regoutputfiles,"AlbedoFinalData.mat");
resdata = cellstr(pwnames);
load(datafname,resdata{:});


% 1. Defined everywhere (except south of 60S)
% ---------------------
% a. Prepare map scale and colorbar    
valcat = [-1200,-800:200:-200,-100,-50,-20,20,50,100,200:200:800,1200];
nbcc = floor((length(valcat)-1)/2);
x = colormap(jet);
lastblue = find(x(:,1)==0,1,'last');
blueint = floor(lastblue/nbcc);
bluescale = x(1:blueint:lastblue,:);
nlb = numel(bluescale)/3;
if nlb > nbcc, bluescale = bluescale([1:nlb-2,nlb],:); end
firstyellow = find(x(:,1)==1,1,'first');
yellowint = floor((length(x)-firstyellow)/nbcc);
yellowscale = x(firstyellow:yellowint:length(x),:);
nly = numel(yellowscale)/3;
if nly > nbcc, yellowscale = yellowscale([1,3:nly],:); end
greyval = [0.7,0.7,0.7];
catblind = cat(1,bluescale,greyval,yellowscale);

% b. Loop through all the land conversions
for pp = 1 : npw
    pname = pwnames(pp);
    ii = find(pwid==pp);
    
    for ll = 1 : length(ii)
        z = ii(ll);
        lc1 = biomes_igbp == pathwayslandcover(z,1);
        lc2 = biomes_igbp == pathwayslandcover(z,2);
        
        % c. Loop through all kernels
        for kk = 1 : nkernels
            eval(strcat("data = ",pname,"(:,:,kk,ll);"))
            data(Antarctica) = NaN;  %#ok<SAGROW>
            val1 = round(min(data,[],'all'));
            val2 = round(max(data,[],'all'));
            
            catdata = discretize(data,valcat);
            
            figure(ll);clf
            h = gcf;
            h.Units = 'pixels';
            h.Position = [1,31,2560,1333];
            axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
            load coast
            pcolorm(lats,lons,catdata)
            caxis([1 length(valcat)])
            colormap(catblind)
            c = colorbar;
            c.Ticks = 1:length(valcat);
            c.TickLabels = num2cell(valcat);
            c.FontSize = 20;
            c.Label.String = 'CO_2e [t CO_2 ha^-^1]';
            c.Ticks = 1:length(valcat);
            c.TickLabels = num2cell(valcat);
            plotm(lat,long,'Color','Black')
            ax = gca;
            ax.Title.FontSize = 24;
            ax.Title.String = {strcat(IGBPBiomes(lc1)," to ",IGBPBiomes(lc2));...
                strcat("\rm CO_2e of landcover change albedo-induced mean radiative forcing @ ToA (",...
                kernels(kk),") - range is [",num2str(round(val1)),",",num2str(round(val2)),"]")};
            gridm
            fname = strcat(resultslocalfolder,"Figures\",pathwayslandcover(z,1),"2",...
                pathwayslandcover(z,2),".",kernels(kk));
            print(strcat(fname,".jpg"),"-djpeg")
        end
    end
    clear pname ii ll z data fname
end
            

% 2. Defined only where initial land cover is present in the 2016 MODIS map
% -------------------------------------------------------------------------
landcovermap = strcat(mapdir,"LandCover\MCD12C1.A2016001.006.2018324172410.hdf");
LCP = hdfread(landcovermap,"Land_Cover_Type_1_Percent");

for pp = 1 : npw
    pname = pwnames(pp);
    ii = find(pwid==pp);
    
    for ll = 1 : length(ii)
        z = ii(ll);
        lc1 = biomes_igbp == pathwayslandcover(z,1);
        lc2 = biomes_igbp == pathwayslandcover(z,2);
        lc1abs = LCP(:,:,lc1) == 0;

        for kk = 1 : nkernels
            eval(strcat("data = ",pname,"(:,:,kk,ll);"))
            
            data(lc1abs) = NaN;     %#ok<SAGROW>    % Eliminate pixels devoid of land cover 1
            val1 = round(min(data,[],'all'));
            val2 = round(max(data,[],'all'));
            
            catdata = discretize(data,valcat);
                        
            figure(10+ll);clf
            h = gcf;
            h.Units = 'pixels';
            h.Position = [1,31,2560,1333];
            axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
            load coast
            pcolorm(lats,lons,catdata)
            caxis([1 length(valcat)])
            colormap(catblind)
            c = colorbar;
            c.Ticks = 1:length(valcat);
            c.TickLabels = num2cell(valcat);
            c.FontSize = 20;
            c.Label.String = 'CO_2e [t CO_2 ha^-^1]';
            c.Ticks = 1:length(valcat);
            c.TickLabels = num2cell(valcat);
            plotm(lat,long,'Color','Black')
            ax = gca;
            ax.Title.FontSize = 24;
            ax.Title.String = {strcat(IGBPBiomes(lc1)," (where present in MODIS) to ",...
                IGBPBiomes(lc2));...
                strcat("\rm CO_2e of landcover change albedo-induced mean radiative forcing @ ToA (",...
                kernels(kk),") - range is [",num2str(round(val1)),",",num2str(round(val2)),"]")};
            gridm
            fname = strcat(resultslocalfolder,"Figures\Present.",pathwayslandcover(z,1),"2",...
                pathwayslandcover(z,2),".",kernels(kk));
            print(strcat(fname,".jpg"),"-djpeg")
        end
    end
    clear pname ii ll z data fname
end




