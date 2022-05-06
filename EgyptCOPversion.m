% NCS - Global Albedo - TNC-Clark Collaboration - Version 1B - Reduced dataset
%
% This code is used to produce a reduced/compacted set of maps from the version 1 (Glasgow COP).
%   It also creates a "most likely forest" maps wich assigns the forest-type with most area in each
%   pixel or when no forest is present, the forest-type that is most prevalent in the neighboring
%   pixels or even the most prevalent in the region (same climate zone, same continent) if none is
%   present in the neighboring pixels. Similarly, it created a "most likely open lands" where open
%   lands are IGBP classes "grassland", "cropland", "open shurbland" and "Cropland/natural 
%   vegetation mosaics".
%
% Created 3/4/22 by Natalia Hasler, Clark University


clearvars



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Define variables and Import data
% ************************************

% 1. Define variables
% -------------------

% ***************
% a. User-defined
dstat = ["median","min","max"];         % Desired statistics on kernels
opportunitybiomes = ...                 % For input current land cover, ignore water, snow&ice,
    ["ENF","EBF","DNF","DBF","MF",...   %   barren and urban as they are unlikely opportunity biomes
    "CSH","OSH","WSA","SAV","GRA","WET","CRO","MOS"];
forestbiomes = ...                      % Forest biomes in output maps (i.e. reforestation pathway)
    ["ENF","EBF","DNF","DBF","MF"];
extforbiomes = ...                      % Forest biomes in input maps (i.e. deforestation pathway)
    ["ENF","EBF","DNF","DBF","MF","WSA","SAV"];
openlandbiomes = ...                    % Open land biomes, used both as input and ouput land covers         
    ["OSH","GRA","CRO","MOS"];
% ***************

% b. Input maps and datasets
landcovermap = "H:\OriginalDatasets\MODISLandCover\MCD12C1.A2020001.006.2021362215328.hdf";
ecosystemsmap = "H:\IntermediaryFiles\WE_005_resample.tif";
welegendfile = "F:\GlobalAlbedo\WorldRegions\WorldTerrestrialEcosystems.xlsx";
regionsmap = "H:\IntermediaryFiles\WR005.tif";
generalalbedodata = "C:\Work\GlobalAlbedo\AlbedoGeneralData.mat";
ecosysandclim = "H:\IntermediaryFiles\EcosystemsAndClimates.mat";
regoutputfiles = "C:\Work\GlobalAlbedo\SubData\";
resultsmatfile = "H:\IntermediaryFiles\EgyptCOP.mat";
resultmaps = "H:\ResultMaps\EgyptCOP\";
figuredir = "C:\Work\GlobalAlbedo\Figures\EgyptCOP\";

% c. Input variables
biomes_igbp = ["WAT","ENF","EBF",...    % IGBP Biomes 3-letter acronyms
    "DNF","DBF","MF","CSH","OSH","WSA","SAV","GRA","WET","CRO","URB","MOS","SNO","BAR"];
IGBPBiomes = ["Water bodies"; ...       % IGBP Biomes names
    "Evergreen needleleaf forests"; "Evergreen broadleaf forests"; "Deciduous needleleaf forests";...
    "Deciduous broadleaf forests"; "Mixed forests"; "Closed shrublands"; "Open shrublands"; ...
    "Woody savannas"; "Savannas"; "Grasslands"; "Permanent wetlands"; "Croplands"; ...
    "Urban and built-up lands"; "Cropland/natural vegetation mosaics"; "Snow and ice"; "Barren"];
regionnames = ["Antarctica","Asiatic Russia","Australia/New Zealand","Caribbean","Central America",...
    "Central Asia","Eastern Africa","Eastern Asia","Eastern Europe","European Russia",...
    "Melanesia","Micronesia","Middle Africa","Northern Africa","Northern America",...
    "Northern Europe","Polynesia","South America","Southeastern Asia","Southern Africa",...
    "Southern Asia","Southern Europe","Western Africa","Western Asia","Western Europe"];
nreg = length(regionnames);


% 2. Import Data
% --------------
% a. Land cover
LCP = hdfread(landcovermap,"Land_Cover_Type_1_Percent");
LC = hdfread(landcovermap,"Majority_Land_Cover_Type_1");

% b. Regions
RE = readgeoraster(regionsmap);

% c. World ecosystems
WE = readgeoraster(ecosystemsmap);
WEinfo = georasterinfo(ecosystemsmap);


% 3. Extract/Define additional information
% ----------------------------------------
% a. Dimentions and scale
[nlat,nlon,nbiomes] = size(LCP);
latlonscale = 180/nlat;
lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[lons,lats] = meshgrid(lon05,lat05);
regnoval = RE(1,1);                     % The info files gives me a strange value
wenoval = WEinfo.MissingDataIndicator;

% b. Define world chunks to accelerate calculations
biggestblock = 20;
blocksize = biggestblock / latlonscale;
nblx = nlon / blocksize;
nbly = nlat / blocksize;
nblocks = nblx * nbly;
presland = true(nblocks,1);             % Presence of land in given world block

% c. No data values
u8noval = 2^8-1;
i16noval = -2^15;

% c. Land Mask (any pixel with at least 1% non-water) & eliminate Antarctica
if LCP(1,1,1) == 0; biomes_igbp = [biomes_igbp(2:length(biomes_igbp)),biomes_igbp(1)]; end
landmask = true(nlat,nlon);
i = biomes_igbp =="WAT";
kk = LCP(:,:,i) == 100;                 % indices of all water pixels
landmask(kk) = false;
noAntarcticalandmask = landmask;
latA = -60;
Antarctica = lats <= latA;
noAntarcticalandmask(Antarctica) = false;
clear kk i latA

% d. Define biomes for different input/output maps
% i. Current land cover
obi = find(ismember(biomes_igbp,opportunitybiomes)) - 1;
% ii. Forests
fp = ismember(biomes_igbp,forestbiomes);
fpi = find(fp) - 1;
nfb = length(forestbiomes);
efpi = find(ismember(biomes_igbp,extforbiomes)) - 1;
forestprop = LCP(:,:,fp);
% iii. Open land biomes
olp = ismember(biomes_igbp,openlandbiomes);
olpi = find(olp)-1;
nolb = length(openlandbiomes);
openlandprop = LCP(:,:,olp);
clear fp olp

% e. Get pixel areas
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


% 4. Match grid size
% ------------------
%    (Don't know why this one did not cover the entire globe ...)
wesize = WEinfo.RasterReference.RasterSize;
welat = round(WEinfo.RasterReference.LatitudeLimits,2);
welon = round(WEinfo.RasterReference.LongitudeLimits,2);
WE05 = NaN(nlat,nlon);
clat05 = 90 : -latlonscale : -90;
clon05 = -180 : latlonscale : 180;
if wesize(1) ~= nlat
    sc = round(WEinfo.RasterReference.CellExtentInLatitude,2);
    if sc ~= latlonscale
        error("world ecosystem scale is not the right one");
    end
    la1 = find(clat05 == max(welat),1);
    la2 = find(clat05 == min(welat),1);
    if (la2 - la1) ~= wesize(1)
        error("Check my world ecosystem grid for latitudes")
    else
        la2 = la2 - 1;
    end    
else
    la1 = 1;
    la2 = nlat;
end
if wesize(2) ~= nlon
    lo1 = find(clon05 == min(welon),1);
    if isempty(lo1)
        lo1 = 1;
        low1 = 1 + (clon05(1) - min(welon))/latlonscale;
    else
        low1 = 1;
    end
    lo2 = find(clon05 == max(welon),1);
    if isempty(lo2)
        lo2 = nlon;
        low2 = wesize(2) - (max(welon) - clon05(length(clon05)))/latlonscale;
    else
        low2 = wesize(2);
        lo2 = lo2-1;
    end
end
 
WE05(la1:la2,lo1:lo2) = WE(:,low1:low2);
WE05(isnan(WE05)) = wenoval;
WE05 = uint16(WE05);


% 5. Extract climate zones from ecosystem file
% --------------------------------------------
opts = detectImportOptions(welegendfile);
vars = string(opts.VariableNames);
T = readtable(welegendfile,opts);
ci = vars == "Temp_Moist";
cllst = string(T{:,ci});
climatelist = unique(cllst);
nclim  = length(climatelist);
climate = ones(nlat,nlon,'uint8') .* u8noval;
weclim = zeros(length(climatelist),length(cllst));
mxl = 0;
for cc = 1 : nclim
    wel = find(strcmp(cllst,climatelist(cc)));
    weclim(cc,1:length(wel)) = wel;
    if length(wel) > mxl; mxl = length(wel);end
    for ww = 1 : length(wel)
        k = wel(ww);
        pos = WE05 == k;
        climate(pos) = cc;
    end
end
weclim = weclim(:,1:mxl);

clear Antarctica REinfo T WE WEinfo ans cc ci cllst computer earthellipsoid ecosystemsmap i k ...
    la1 la2 landcovermap lat1 lat2 latA lats lo1 lo2 lon1 lon2 low1 low2 mxl opts pix pos ...
    regionsmap sc wel welat welegendfile welon wesize ww

save(ecosysandclim,'RE','WE05','climate','climatelist','regionnames','weclim');


% 6. Define output variables
% --------------------------
% a. Input land cover maps and tables
inmap = ones(nlat,nlon,'uint8') .* u8noval;
currentlandcover = inmap;               % All non-water land cover (only for reference)
currentcleanlandcover = inmap;          % Current land cover, excl. water, urban, snow&ice, barren
currentopenland = inmap;                % Current open lands map
currentforest = inmap;                  % Current forest map (including savannas)
LCPix = zeros(nbiomes+1,4);             % Comparative table for those four maps

% b. Output land cover maps and tables, i.e. "Most likely forest/open land"
RegionalMLF = zeros(nreg,nclim);        % Table of most likely forest given region/climate
RegionalMLOL = zeros(nreg,nclim);       % Table of most likely open land class given region/climate
MLOL = inmap;                           % Map of most likely open land class
MLF = inmap;                            % Map of most likely
MLC = NaN(180,nreg+1);                  % Most likely climate per latitude/region
pixregclim = zeros(nreg,nclim);         % Table of pixel count in each region/climate
arearegclim = zeros(nreg,nclim);        % Table of global area in each region/climate

% b. Albedo effects maps
in = ["CLC","FOR","OPL"];
out = ["FOR","OPL","URB"];
ine = ["CurrentLandCover","CurrentForest","CurrentOpenLands"];
inl = ["Current Land Cover","Current Forest","Current Open Lands"];
outl = ["Forest","Open Lands","Urban"];
oute = ["Forest","OpenLand","Urban"];
if sum(strcmp(dstat,["median","min","max"])) == 3
    st = ["med","min","max"];
    ste = ["median","minimum","maximum"];
else
    error("correct map names")
end
outmap = ones(nlat,nlon,'int16') .* i16noval; %#ok<NASGU>
outvarnames = strings((3*3-2)*3,1);
outmapnames = strings((3*3-2)*3,1);
maplegends = strings((3*3-2)*3,2);
ct = 0;
for ii = 1 : length(in)
    for oo = 1 : length(out)
        if strcmp(in(ii),out(oo)) == 0
            for ss = 1 : length(st)
                ct = ct + 1;
                vn = strcat(in(ii),out(oo),st(ss)); %#ok<NASGU>
                leg = strcat(inl(ii)," to ",outl(oo));
                leg2 = ste(ss);
%                 eval(strcat("vn = '",in(ii),out(oo),st(ss),"';"));
                eval(strcat(in(ii),out(oo),st(ss)," = outmap;"));              
                eval("outvarnames(ct) = string(vn);");
                eval(strcat("outmapnames(ct) = string('",ine(ii),"2",oute(oo),"_",ste(ss),"');"))
                eval("maplegends(ct,1) = string(leg);")
                eval("maplegends(ct,2) = string(leg2);")
            end
        end
    end 
end
clear inmap outmap in out st vn ii oo ss ct ine oute ste inl outl

save(resultsmatfile,'IGBPBiomes','biggestblock','biomes_igbp','blocksize','efpi','extforbiomes',...
    'figuredir','forestbiomes','fpi','globalpixelarea','i16noval','landmask','lat05','latlonscale',...
    'lon05','lons','nbiomes','nblocks','nblx','nbly','nclim','climatelist','nfb','nlat','nlon',...
    'noAntarcticalandmask','nolb','nreg','obi','olpi','openlandbiomes','opportunitybiomes',...
    'outvarnames','outmapnames','maplegends','regionnames','u8noval','regoutputfiles','resultmaps')




%% B. Create initial land cover maps
% **********************************

% 1. Gap-fill climate zones and regions
% -------------------------------------
% a. regions
[misla,mislo] = find(((RE == regnoval) + noAntarcticalandmask) == 2);
for k = 1 : numel(misla)
    RE = geogapfill(misla(k),mislo(k),RE,regionnames,regnoval,globalpixelarea);
end
if regnoval ~= u8noval
    RE(RE==regnoval) = u8noval;             % for consistency between datasets
end
clear misla mislo k

% b. Climate
% i. Create a most likely climate per latitude (and continent/region, if applicable)
for ll = 1 : 180
    i = (ll-1)*1/latlonscale + 1 : ll*1/latlonscale;
    clim = double(climate(i,:));
    clim(clim==u8noval) = NaN;
    reg = RE(i,:);
    area = globalpixelarea(i,:);
    if sum(isfinite(clim),'all')>0
        for rr = 1 : nreg + 1
            if rr == nreg + 1
                ii = true(1/latlonscale,nlon);
            else
                ii = reg == rr;
            end
            if sum(isfinite(clim(ii)),'all') > 0
                s = min(clim(ii),[],'all','omitnan'):max(clim(ii),[],'all','omitnan');
                n = histcounts(clim(ii));
                if sum(n == max(n)) == 1
                    val = s(n==max(n));
                    MLC(ll,rr) = val;
                else
                    A = zeros(length(s),1);
                    smarea = area(ii);
                    for k = 1 : length(s)
                        pp = clim(ii) == s(k);
                        A(k) = sum(smarea(pp));
                    end
                    kk = find(A == max(A));
                    if numel(kk) == 1
                        MLC(ll,rr) = s(kk);
                    end
                end
            end
        end
    end
end
MLC(isnan(MLC)) = u8noval;
MLC = uint8(MLC);
clear ll i clim reg area rr ii s n val A smarea k pp kk

% ii. Gap fill
[misla,mislo] = find(((climate == u8noval) + noAntarcticalandmask) == 2);
for k = 1 : numel(misla)
    climate = geogapfill(misla(k),mislo(k),climate,climatelist,u8noval,globalpixelarea,MLC,RE);
end
clear misla mislo k

save(ecosysandclim,'RE','climate','-append');
save(resultsmatfile,'MLC','-append');


% 2. Create current land (for reference only)
% ----------------------
%    Current land cover is the MODIS majority land cover layer with added values on partial water
%    pixels (i.e. on all landmask pixels, defined as >1% non-water)
blcurlc = ones(blocksize,blocksize,nblocks) .* u8noval;
for bb = 1 : nblocks
    y = ceil(bb/nblx);              % subregion index in the y direction
    x = bb - (y-1)*nblx;            % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    
    ldm = noAntarcticalandmask(i,j);
    if sum(ldm,'all') == 0
        presland(bb) = false;
        continue
    end

    smbl = LC(i,j);
    smbl(ldm==0) = u8noval;
    [a,b] = find(((ldm==1) + (smbl==0)) == 2);
    
    if numel(a) > 0
        blc = bestlandcover(a,b,LCP(i,j,2:17),(1:16),climate(i,j),RE(i,j),globalpixelarea(i,j));
        px = sub2ind([blocksize,blocksize],a,b);
        smbl(px) = blc;
    end
    smbl(smbl==0) = u8noval;
    blcurlc(:,:,bb) = uint8(smbl);
    currentlandcover(i,j) = uint8(smbl);
end
LCPix(:,1) = histcounts(blcurlc,[0:nbiomes,256]);


% 3. Create current "cleaned" land cover (where snow%ice, bare, and urban are removed)
% --------------------------------------
blclcrlc = ones(blocksize,blocksize,nblocks,'uint8') .* u8noval;
for bb = 1 : nblocks
    if presland(bb) == 0, continue; end
    y = ceil(bb/nblx);              % subregion index in the y direction
    x = bb - (y-1)*nblx;            % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    smbl = ones(blocksize,blocksize,'uint8') .* u8noval;
    ldm = noAntarcticalandmask(i,j);
    [a,b] = find(ldm);
    blc = bestlandcover(a,b,LCP(i,j,obi+1),obi,climate(i,j),RE(i,j),globalpixelarea(i,j));
    smbl(ldm) = blc;
    blclcrlc(:,:,bb) = uint8(smbl);
    currentcleanlandcover(i,j) = uint8(smbl);
end
LCPix(:,2) = histcounts(blclcrlc,[0:nbiomes,256]);

% c. Initial as majority of open land
%    (i.e. open shrubland (7), grassland (10), croplands (12), cropland/natural veg mosaics (14))
blopenlc = ones(blocksize,blocksize,nblocks) .* u8noval;
for bb = 1 : nblocks
    if presland(bb) == 0, continue; end
    y = ceil(bb/nblx);              % subregion index in the y direction
    x = bb - (y-1)*nblx;            % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    smbl = zeros(blocksize,blocksize);
    ldm = noAntarcticalandmask(i,j);
    [a,b] = find(ldm);
    blc = bestlandcover(a,b,LCP(i,j,olpi+1),olpi,climate(i,j),RE(i,j),globalpixelarea(i,j));
    smbl(ldm==1) = blc;
    smbl(ldm==0) = u8noval;
    blopenlc(:,:,bb) = uint8(smbl);
    currentopenland(i,j) = uint8(smbl);
end
LCPix(:,3) = histcounts(blopenlc,[0:nbiomes,256]);

% c. Initial as majority of forest (including savannas)
blforlc = ones(blocksize,blocksize,nblocks) .* u8noval;
for bb = 1 : nblocks
    if presland(bb) == 0, continue; end
    y = ceil(bb/nblx);              % subregion index in the y direction
    x = bb - (y-1)*nblx;            % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    smbl = zeros(blocksize,blocksize);
    ldm = noAntarcticalandmask(i,j);
    [a,b] = find(ldm);
    blc = bestlandcover(a,b,LCP(i,j,efpi+1),efpi,climate(i,j),RE(i,j),globalpixelarea(i,j));
    smbl(ldm==1) = blc;
    smbl(ldm==0) = u8noval;
    blforlc(:,:,bb) = uint8(smbl);
    currentforest(i,j) = uint8(smbl);
end
LCPix(:,4) = histcounts(blforlc,[0:nbiomes,256]);
clear bb x y i j smbl ldm a b blc nsmbl

save(resultsmatfile,'LCPix','blcurlc','blclcrlc','blopenlc','blforlc','currentforest',...
    'currentopenland','currentcleanlandcover','currentlandcover','presland','-append')




%% C. Create Most likely forest/open lands maps
% *********************************************

% 1. Create a table of most likely forest/open lands per region/climate zone
% --------------------------------------------------------------------------
%    (When several biomes have the same likelyhood, create an index with powers of 10)
% a. Loop through climate/regions to get most prevalent biome
index = [10,100,1000,10000,100000];
for rr = 1 : nreg
    sr = RE == rr;
    for cc = 1 : nclim
        A = zeros(nfb,1);
        A2 = zeros(nolb,1);
        sc = climate == cc;
        ss = sr+sc == 2;
        if sum(ss,'all') > 0
            pixregclim(rr,cc) = sum(ss,'all');
            arearegclim(rr,cc) = sum(globalpixelarea(ss));
            for ff = 1 : nfb
                F = forestprop(:,:,ff);
                A(ff) = sum(double(F(ss)) .* globalpixelarea(ss));
            end            
            if max(A) > 0 && sum(A == max(A)) == 1
                RegionalMLF(rr,cc) = find(A==max(A),1);
            elseif max(A) > 0
                RegionalMLF(rr,cc) = sum(index(A==max(A)));
            end
            for oo = 1 : nolb
                O = openlandprop(:,:,oo);
                A2(oo) = sum(double(O(ss)) .* globalpixelarea(ss));
            end            
            if max(A2) > 0 && sum(A2 == max(A2)) == 1
                RegionalMLOL(rr,cc) = find(A2==max(A2),1);
            elseif max(A2) > 0
                RegionalMLOL(rr,cc) = sum(index(A2==max(A2)));
            end
            
        end
    end
end
clear index rr sr cc A A2 sc ss ff F oo O

% b. Adjust open land index to reflect biome category
for oo = 1 : nolb
    RegionalMLOL(RegionalMLOL == oo) = olpi(oo);
end

save(ecosysandclim,'pixregclim','arearegclim','-append');
save(resultsmatfile,'RegionalMLF','RegionalMLOL','-append');


% 2. Create several resolution gridded maps for both forests and open lands
% -------------------------------------------------------------------------
PF01  = ones(nlat,nlon,'uint8') *u8noval;
PF025 = ones(nlat,nlon,'uint8') *u8noval;
PF05  = ones(nlat,nlon,'uint8') *u8noval;
PF1   = ones(nlat,nlon,'uint8') *u8noval;
PF25  = ones(nlat,nlon,'uint8') *u8noval;
PF5   = ones(nlat,nlon,'uint8') *u8noval;
PF10  = ones(nlat,nlon,'uint8') *u8noval;
PF20  = ones(nlat,nlon,'uint8') *u8noval;
PFblock20 = zeros(nreg,nclim,nblx,nbly);

POL01  = ones(nlat,nlon,'uint8') *u8noval;
POL025 = ones(nlat,nlon,'uint8') *u8noval;
POL05  = ones(nlat,nlon,'uint8') *u8noval;
POL1   = ones(nlat,nlon,'uint8') *u8noval;
POL25  = ones(nlat,nlon,'uint8') *u8noval;
POL5   = ones(nlat,nlon,'uint8') *u8noval;
POL10  = ones(nlat,nlon,'uint8') *u8noval;
POL20  = ones(nlat,nlon,'uint8') *u8noval;
POLblock20 = zeros(nreg,nclim,nblx,nbly);

% b. Loop through blocks to calculate prevalences at different scales
for bb = 1 : nblocks
    y = ceil(bb/nblx);              % subregion index in the y direction
    x = bb - (y-1)*nblx;            % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    
    zonelandmask = noAntarcticalandmask(i,j);
    if sum(zonelandmask,'all') == 0
        presland(bb) = false;
        continue
    end
    zonefor = forestprop(i,j,:);
    zoneopen = openlandprop(i,j,:);
    zonereg = RE(i,j);
    zoneclim = climate(i,j);
    area = globalpixelarea(i,j);

    
    % Calculate prevalence over .1 deg resolution grid
    resolution = 0.1;
    [BPF01,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL01,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over .25 deg resolution grid
    resolution = 0.25;
    [BPF025,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL025,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over .5 deg resolution grid
    resolution = 0.5;
    [BPF05,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL05,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over 1 deg resolution grid
    resolution = 1;
    [BPF1,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL1,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over 2.5 deg resolution grid
    resolution = 2.5;
    [BPF25,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL25,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over 5*5 deg resolution grid
    resolution = 5;
    [BPF5,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL5,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over 10 deg resolution grid
    resolution = 10;
    [BPF10,~,~,~] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL10,~,~,~] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Calculate prevalence over entire block
    resolution = 20;
    [BPF20,ftbl,frlst,fclst] = prevalentbiomegrid(resolution,zonefor,zonereg,zoneclim,area,RegionalMLF,...
        zonelandmask,fpi,latlonscale,u8noval);
    [BPOL20,olbl,olrlst,olclst] = prevalentbiomegrid(resolution,zoneopen,zonereg,zoneclim,area,RegionalMLOL,...
        zonelandmask,olpi,latlonscale,u8noval);
    
    % Save values
    PF01(i,j) = BPF01;
    PF025(i,j) = BPF025;
    PF05(i,j) = BPF05;
    PF1(i,j) = BPF1;
    PF25(i,j) = BPF25;
    PF5(i,j) = BPF5;
    PF10(i,j) = BPF10;
    PF20(i,j) = BPF20;
    
    PFblock20(frlst,fclst,x,y) = ftbl;
    
    POL01(i,j) = BPOL01;
    POL025(i,j) = BPOL025;
    POL05(i,j) = BPOL05;
    POL1(i,j) = BPOL1;
    POL25(i,j) = BPOL25;
    POL5(i,j) = BPOL5;
    POL10(i,j) = BPOL10;
    POL20(i,j) = BPOL20;
    
    POLblock20(olrlst,olclst,x,y) = olbl;
    
end

save(resultsmatfile,'presland','PF01','PF025','PF05','PF1','PF25','PF5','PF10','PF20','POL01',...
    'POL025','POL05','POL1','POL25','POL5','POL10','POL20','PFblock20','POLblock20','-append');
    

% 3. Populate map where forest/openland  is present (it is way faster if I chunk it up ...)
% -------------------------------------------------
for bb = 1 : nblocks
    y = ceil(bb/nblx);                  % subregion index in the y direction
    x = bb - (y-1)*nblx;                % subregion index in the x derection
    j = (x-1)*blocksize + 1 : x*blocksize;
    i = (y-1)*blocksize + 1 : y*blocksize;
    smbl = zeros(blocksize,blocksize);
    ldm = noAntarcticalandmask(i,j);
    [a,b] = find(ldm);
    blc = bestlandcover(a,b,forestprop(i,j,:),fpi,climate(i,j),RE(i,j),globalpixelarea(i,j));
    smbl(ldm==1) = blc;
    smbl(ldm==0) = u8noval;
    MLF(i,j) = uint8(smbl);
    smbl = zeros(blocksize,blocksize);
    blc = bestlandcover(a,b,openlandprop(i,j,:),olpi,climate(i,j),RE(i,j),globalpixelarea(i,j));
    smbl(ldm==1) = blc;
    smbl(ldm==0) = u8noval;
    MLOL(i,j) = uint8(smbl);   
end


% 4. Populate the rest of the grid
% --------------------------------
%    To do so, I loop through the different resolution grids until I find a value.

% a. Fill forest cover
[misla,mislo] = find(MLF == 0);

for px = 1 : numel(misla)
    la = misla(px);
    lo = mislo(px);
    val = PF01(la,lo);
    if val == 0
        val = PF025(la,lo);
        if val == 0
            val = PF05(la,lo);
            if val == 0
                val = PF1(la,lo);
                if val == 0
                    val = PF25(la,lo);
                    if val == 0
                        val = PF5(la,lo);
                        if val == 0
                            val = PF10(la,lo);
                            if val == 0
                                val = PF20(la,lo);
                                if val == 0
                                    x = ceil(lo/blocksize);
                                    y = ceil(la/blocksize);
                                    a = x-1:x+1; a(a<1) = NaN; a(a>nblx) = NaN; a = a(isfinite(a));
                                    b = y-1:y+1; b(b<1) = NaN; b(b>nbly) = NaN; b = b(isfinite(b));
                                    clim = climate(la,lo);
                                    reg = RE(la,lo);
                                    spr = PFblock20(reg,clim,a,b);
                                    n = histcounts(spr,1:6);
                                    if max(n) == 0
                                        val = RegionalMLF(reg,clim);
                                    elseif sum(n==max(n)) == 1
                                        val = find(n==max(n),1);
                                    else
                                        f2 = find(n==max(n));
                                        if ismember(RegionalMLF(reg,clim),f2)
                                            val = RegionalMLF(reg,clim);
                                        else
                                            val = f2(random('Discrete Uniform',numel(f2)));
                                        end
                                    end
                                elseif isnan(val)
                                    error("I have discrepancy between my 20res grid and forest landmask")
                                end
                            elseif isnan(val)
                                error("I have discrepancy between my 10res grid and forest landmask")
                            end
                        elseif isnan(val)
                            error("I have discrepancy between my 5res grid and forest landmask")
                        end
                    elseif isnan(val)
                        error("I have discrepancy between my 2.5res grid and forest landmask")
                    end
                elseif isnan(val)
                    error("I have discrepancy between my 1res grid and forest landmask")
                end
            elseif isnan(val)
                error("I have discrepancy between my 0.5res grid and forest landmask")
            end
        elseif isnan(val)
            error("I have discrepancy between my 0.25res grid and forest landmask")
        end
    elseif isnan(val)
        error("I have discrepancy between my 0.1res grid and forest landmask")
    end
    MLF(la,lo) = val;
end

% b. Fill open land cover
[misla,mislo] = find(MLOL == 0);

for px = 1 : numel(misla)
    la = misla(px);
    lo = mislo(px);
    val = POL01(la,lo);
    if val == 0
        val = POL025(la,lo);
        if val == 0
            val = POL05(la,lo);
            if val == 0
                val = POL1(la,lo);
                if val == 0
                    val = POL25(la,lo);
                    if val == 0
                        val = POL5(la,lo);
                        if val == 0
                            val = POL10(la,lo);
                            if val == 0
                                val = POL20(la,lo);
                                if val == 0
                                    x = ceil(lo/blocksize);
                                    y = ceil(la/blocksize);
                                    a = x-1:x+1; a(a<1) = NaN; a(a>nblx) = NaN; a = a(isfinite(a));
                                    b = y-1:y+1; b(b<1) = NaN; b(b>nbly) = NaN; b = b(isfinite(b));
                                    clim = climate(la,lo);
                                    reg = RE(la,lo);
                                    spr = POLblock20(reg,clim,a,b);
                                    n = histcounts(spr,1:6);
                                    if max(n) == 0
                                        val = RegionalMLOL(reg,clim);
                                    elseif sum(n==max(n)) == 1
                                        val = find(n==max(n),1);
                                    else
                                        f2 = find(n==max(n));
                                        if ismember(RegionalMLOL(reg,clim),f2)
                                            val = RegionalMLOL(reg,clim);
                                        else
                                            val = f2(random('Discrete Uniform',numel(f2)));
                                        end
                                    end
                                elseif isnan(val)
                                    error("I have discrepancy between my 20res grid and forest landmask")
                                end
                            elseif isnan(val)
                                error("I have discrepancy between my 10res grid and forest landmask")
                            end
                        elseif isnan(val)
                            error("I have discrepancy between my 5res grid and forest landmask")
                        end
                    elseif isnan(val)
                        error("I have discrepancy between my 2.5res grid and forest landmask")
                    end
                elseif isnan(val)
                    error("I have discrepancy between my 1res grid and forest landmask")
                end
            elseif isnan(val)
                error("I have discrepancy between my 0.5res grid and forest landmask")
            end
        elseif isnan(val)
            error("I have discrepancy between my 0.25res grid and forest landmask")
        end
    elseif isnan(val)
        error("I have discrepancy between my 0.1res grid and forest landmask")
    end
    MLOL(la,lo) = val;
end

save(resultsmatfile,'MLF','MLOL','-append');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% C. Create pathway maps
%  **********************

% 1. Define global variables
% --------------------------
load(generalalbedodata,'nmonths','nkernels','kernels','kernelscale','albsn')
monthln = ["January","February","March","April","May","June","July","August",...
    "September","October","November","December"];
kervars = cellstr(lower(kernels));
albvar = cellstr(albsn);
othvars = cellstr(["snowcover","diffusefraction","pixelarea","regi","regj"]);
emptymap = ones(blocksize,blocksize,'int16') .* i16noval;


% 2. Loop through world "bigblocks" (20x20deg chunks)
% ---------------------------------
for bb = 1 : nblocks
    tic
    if presland(bb) == true
        
        
        % 3. Get/define local variables
        % -----------------------------
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(bb),".mat");
        load(subdatafname,kervars{:},albvar{:},othvars{:});
        forest = MLF(regi,regj);
        openland = MLOL(regi,regj);
        BOL = blopenlc(:,:,bb);
        BEF = blforlc(:,:,bb);
        BLC = blclcrlc(:,:,bb);            
        blocklandmask = noAntarcticalandmask(regi,regj);        
        
        clcformed = emptymap; clcoplmed = emptymap; clcurbmed = emptymap;
        clcformin = emptymap; clcoplmin = emptymap; clcurbmin = emptymap;
        clcformax = emptymap; clcoplmax = emptymap; clcurbmax = emptymap;
        
        foroplmed = emptymap; forurbmed = emptymap;
        foroplmin = emptymap; forurbmin = emptymap;
        foroplmax = emptymap; forurbmax = emptymap;
        
        oplformed = emptymap; oplurbmed = emptymap;
        oplformin = emptymap; oplurbmin = emptymap;
        oplformax = emptymap; oplurbmax = emptymap;        
        
        ker = NaN(blocksize,blocksize,12,nkernels);
        for kk = 1 : nkernels
            eval(strcat("ker(:,:,:,kk) = ",lower(kernels(kk))," .* kernelscale(kk);"))
        end
        
        
        % 4. Loop through land pixels
        % ---------------------------
        [a,b] = find(blocklandmask);
        
        for px = 1 : numel(a)
            i = a(px);
            j = b(px);
            area = pixelarea(i,j);
            nsbs = squeeze(nosnowblacksky(i,j,:,:));
            nsws = squeeze(nosnowwhitesky(i,j,:,:));
            scbs = squeeze(snowcovblacksky(i,j,:,:));
            scws = squeeze(snowcovwhitesky(i,j,:,:));
            snow = squeeze(snowcover(i,j,:));
            diffuse = squeeze(diffusefraction(i,j,:));
            kerval = squeeze(ker(i,j,:,:));
            
            % a. Current land cover maps
            lc = BLC(i,j);
            if ismember(lc,obi)
                lc1 = lc + 1;
                
                % (1) Current land cover to forest (exclude initial five forest classes)
                mlfrst = forest(i,j);
                if ismember(mlfrst,fpi)
                    lc2 = mlfrst + 1;
                    if ismember(lc,fpi) == 0
                        pval = ...
                            ptradforcing(lc1,lc2,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);                        
                        clcformed(i,j) = pval(1);
                        if pval(2) < pval(3)  % Because I change the sign, I might have to swap min/max
                            clcformin(i,j) = pval(2);
                            clcformax(i,j) = pval(3);
                        else
                            clcformin(i,j) = pval(3);
                            clcformax(i,j) = pval(2);
                        end
                        clear pval(1) pval(2) pval(3)
                    end
                end
                
                % (2) Current land cover to open land (exclude initial open land)
                mlopld = openland(i,j);
                if ismember(mlopld,olpi)
                    lc3 = mlopld + 1;
                    if ismember(lc,olpi) == 0
                        pval = ...
                            ptradforcing(lc1,lc3,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                        clcoplmed(i,j) = pval(1);
                        if pval(2) < pval(3)
                            clcoplmin(i,j) = pval(2);
                            clcoplmax(i,j) = pval(3);
                        else
                            clcoplmin(i,j) = pval(3);
                            clcoplmax(i,j) = pval(2);
                        end
                        clear pval(1) pval(2) pval(3)
                    end
                end
                
                % (3) Current land cover to urban
                lc4 = find(biomes_igbp == "URB",1);
                if lc1 ~= lc4
                    pval = ...
                        ptradforcing(lc1,lc4,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                    clcurbmed(i,j) = pval(1);
                    if pval(2) < pval(3)
                        clcurbmin(i,j) = pval(2);
                        clcurbmax(i,j) = pval(3);
                    else
                        clcurbmin(i,j) = pval(3);
                        clcurbmax(i,j) = pval(2);
                    end
                    clear pval(1) pval(2) pval(3)
                end
            end
            
            % b. Deforestation scenarios (from forests and savannas)
            frst = BEF(i,j);
            if frst == lc
                foroplmed(i,j) = clcoplmed(i,j); forurbmed(i,j) = clcurbmed(i,j);
                foroplmin(i,j) = clcoplmin(i,j); forurbmin(i,j) = clcurbmin(i,j);
                foroplmax(i,j) = clcoplmax(i,j); forurbmax(i,j) = clcurbmax(i,j);
                
            elseif ismember(frst,efpi)
                lc1 = frst + 1;
                
                % (1) Forest to open land
                pval = ...
                    ptradforcing(lc1,lc3,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                foroplmed(i,j) = pval(1);
                if pval(2) < pval(3)
                    foroplmin(i,j) = pval(2);
                    foroplmax(i,j) = pval(3);
                else
                    foroplmin(i,j) = pval(3);
                    foroplmax(i,j) = pval(2);
                end
                clear pval(1) pval(2) pval(3)
                                
                % (2) Forest to urban
                pval = ...
                    ptradforcing(lc1,lc4,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                forurbmed(i,j) = pval(1);
                if pval(2) < pval(3)
                    forurbmin(i,j) = pval(2);
                    forurbmax(i,j) = pval(3);
                else
                    forurbmin(i,j) = pval(3);
                    forurbmax(i,j) = pval(2);
                end
                clear pval(1) pval(2) pval(3)
                            
            end
            
            % c. Reforestation scenario (from open lands) and conversion from open land to urban
            opld = BOL(i,j);
            if opld == lc
                oplformed(i,j) = clcformed(i,j); oplurbmed(i,j) = clcurbmed(i,j);
                oplformin(i,j) = clcformin(i,j); oplurbmin(i,j) = clcurbmin(i,j);
                oplformax(i,j) = clcformax(i,j); oplurbmax(i,j) = clcurbmax(i,j);
                
            elseif ismember(opld,olpi)
                lc1 = opld + 1;
                
                % (1) Open land to forest
                pval = ...
                    ptradforcing(lc1,lc2,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                oplformed(i,j) = pval(1);
                if pval(2) < pval(3)
                    oplformin(i,j) = pval(2);
                    oplformax(i,j) = pval(3);
                else
                    oplformin(i,j) = pval(3);
                    oplformax(i,j) = pval(2);
                end
                clear pval(1) pval(2) pval(3)
                
                % (2) Open land to urban
                pval = ...
                    ptradforcing(lc1,lc4,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,dstat);
                oplurbmed(i,j) = pval(1);
                if pval(2) < pval(3)
                    oplurbmin(i,j) = pval(2);
                    oplurbmax(i,j) = pval(3);
                else
                    oplurbmin(i,j) = pval(3);
                    oplurbmax(i,j) = pval(2);
                end
                clear pval(1) pval(2) pval(3)

            end
                        
        end
        
        
        % 5. Save block maps into world maps
        % ----------------------------------
        for vv = 1 : length(outvarnames)
            vname = outvarnames(vv);
            v2 = lower(outvarnames(vv));
            eval(strcat(vname,"(regi,regj) = ",v2,";"));
        end
    
        clear clcformed clcoplmed clcurbmed clcformin clcoplmin clcurbmin clcformax clcoplmax ...
            clcurbmax foroplmed forurbmed foroplmin forurbmin foroplmax forurbmax ...
            oplformed oplurbmed oplformin oplurbmin oplformax oplurbmax frst opld mlfrst mlopld ...
            lc lc1 lc2 lc3 lc4 a b px i j area nsbs nsws scbs scws snow diffuse ...
            kermedval kerminval kermaxval forest openland BLC BOL BEF blocklandmask ...
            subdatafname medker minker maxker
        clear(kervars{:},albvar{:},othvars{:})
    
    end
    toc
    strcat("Done with subregion ",num2str(bb))

end

clear albsn albvar bb emptymap kervars othvars emptymap
savevars = cellstr(outvarnames);
save(resultsmatfile,savevars{:},'-append');

        


%% D. Save maps in GeoTIFF format
%  ******************************

% 1. Define geotiff
% -----------------
latlim = [-90,90];
lonlim = [-180,180];
rastersize = [nlat,nlon];
R = georefcells(latlim,lonlim,rastersize,'ColumnsStartFrom','north');
cmptag.Compression = 'LZW';             % LZW compression is better than the default ...


% 2. Export input data
% --------------------
% a. Climate and regions
fname = strcat(resultmaps,"WorldRegions.tif");
geotiffwrite(fname,RE,R,'TiffTags',cmptag);
fname = strcat(resultmaps,"WorldClimates.tif");
geotiffwrite(fname,climate,R,'TiffTags',cmptag);

% b. Initial land cover
fname = strcat(resultmaps,"CurrentLandCover.tif");
geotiffwrite(fname,currentcleanlandcover,R,'TiffTags',cmptag);
fname = strcat(resultmaps,"CurrentForest.tif");
geotiffwrite(fname,currentforest,R,'TiffTags',cmptag);
fname = strcat(resultmaps,"CurrentOpenLand.tif");
geotiffwrite(fname,currentopenland,R,'TiffTags',cmptag);

% c. final land cover
fname = strcat(resultmaps,"MostLikelyForest.tif");
geotiffwrite(fname,MLF,R,'TiffTags',cmptag);
fname = strcat(resultmaps,"MostLikelyOpenLand.tif");
geotiffwrite(fname,MLOL,R,'TiffTags',cmptag);


% 3. Export output data
% ---------------------
for vv = 1 : length(outvarnames)
    vname = outvarnames(vv);
    eval(strcat("vval = ",vname,";"))
    fname = strcat(resultmaps,outmapnames(vv),".tif");
    geotiffwrite(fname,vval,R,'TiffTags',cmptag);
end





%% E. Export data as tables and figures
% *************************************

% 1. Save tables in excel
% ----------------------
% a. Save current land cover # of pixels
catnames = ["Blank";IGBPBiomes(2:17);"Ocean/Water/Antarctica"];
CurrentLandCoverPix = LCPix(:,1);
ReducedLandCoverPix = LCPix(:,2);
CurrentForestPix = LCPix(:,4);
OpenLandsPix = LCPix(:,3);
T = table(catnames,CurrentLandCoverPix,ReducedLandCoverPix,CurrentForestPix,OpenLandsPix);
writetable(T,strcat(resultmaps,"PixelCounts.xlsx"));
clear T

% b. Save region/climate areas
climorder = [15:-1:10,18,17,16,6:-1:1,9,8,7];
climnames = climatelist(climorder);
T = table('Size',[nclim,nreg+1],'VariableTypes',cat(1,"string",repmat("double",[nreg,1])));
T(:,1) = cellstr(climnames);
T(:,2:nreg+1) = num2cell(permute(arearegclim(:,climorder),[2,1])/10^6);
vnames = cat(2,"Climate",replace(replace(regionnames," ",""),"/",""));
T.Properties.VariableNames = cellstr(vnames);
writetable(T,strcat(resultmaps,"ClimateRegionArea.xlsx"));
clear T

% c. Save most likely forest
regorder = [15,5,4,18,16,25,22,9,10,14,23,13,20,7,24,21,6,2,8,19,3,11,12,17,1];
T = table('Size',[nreg,nclim+1],'VariableTypes',cat(1,"string",repmat("double",[nclim,1])));
T(:,1) = cellstr(permute(regionnames(regorder),[2,1]));
T(:,2:nclim+1) = num2cell(RegionalMLF(regorder,climorder));
vnames = cat(2,"Region",permute(replace(climnames," ",""),[2,1]));
T.Properties.VariableNames = cellstr(vnames);
writetable(T,strcat(resultmaps,"MostLikelyForest.xlsx"));
clear T

% c. Save most likely open land
T = table('Size',[nreg,nclim+1],'VariableTypes',cat(1,"string",repmat("double",[nclim,1])));
T(:,1) = cellstr(permute(regionnames(regorder),[2,1]));
T(:,2:nclim+1) = num2cell(RegionalMLOL(regorder,climorder));
vnames = cat(2,"Region",permute(replace(climnames," ",""),[2,1]));
T.Properties.VariableNames = cellstr(vnames);
writetable(T,strcat(resultmaps,"MostLikelyOpenLand.xlsx"));
clear T


% 2. Print base-layer maps as figures
% -----------------------------------
load coastlines
hotmap = colormap(hot);
jetmap = colormap(jet);
lastblue = find(jetmap(:,1)==0,1,'last');
firstyellow = find(jetmap(:,1)==1,1,'first');

colorblindpalette = [0,0,0;0,73,73;0,146,146;255,109,182;255,182,219;73,0,146;0,109,219;...
    182,109,255;109,182,255;182,219,255;146,0,0;146,73,0;219,109,0;36,255,36;255,255,219;...
    255,255,109]/255;
reordcbp = colorblindpalette([1:3,9:10,6:7,4:5,14,11,8,12:13,15:16],:);
lightgreyval = [.9,.9,.9];
greyval = [0.7,0.7,0.7];

% (1) regions
figure(1); clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
reg = double(RE); reg(reg==255) = NaN;
reg(reg==4) = 26;
reg(reg==10) = 4;
pcolorm(lats,lons,reg)
plotm(coastlat,coastlon,'Color','Black')
colormap(reordcbp)
c = colorbar;
c.Ticks = 1+25/32:25/16:26;
rn = permute(regionnames,[2,1]);
regn = cat(1,strcat(rn(1)," & ",rn(2)),strcat(rn(3)," & Eur. Russia"),rn(5),strcat(rn(6)," & ",rn(7)),...
    rn([8,9,11]),strcat(rn(13)," & ",rn(12)),strcat(rn(15)," & ",rn(14)),rn(16),...
    strcat(rn(17)," & ",rn(18)),rn(19),strcat(rn(20)," & ",rn(21)),...
    rn(22),strcat(rn(23)," & ",rn(24)),strcat(rn(4)," & ",rn(25)));
c.TickLabels = cellstr(regn);
c.FontSize = 18;
ax = gca;
ax.Position = [0.0156,0.0156,0.7511,0.899];
ax.Title.FontSize = 24;
ax.Title.String = "Regions";
gridm
fname = strcat(figuredir,"regions.jpg");
print(fname,"-djpeg")

% (2) climates
clim = double(climate);
clim(clim==u8noval) = NaN;
figure(2);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
% h.Position = [1,31,1920,1093];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
load coastlines
pcolorm(lats,lons,clim)
plotm(coastlat,coastlon,'Color','Black')
clcm = [colorblindpalette(1:5,:);1,219/255,1;colorblindpalette(6:10,:);219/255,0,0;colorblindpalette(11:16,:)];
colormap(clcm)
c = colorbar;
c.Ticks = 1+17/36:17/18:18;
c.TickLabels = cellstr(climatelist);
c.FontSize = 18;
ax = gca;
ax.Position = [0.005,0.1725,1,0.7125];
ax.Title.FontSize = 24;
ax.Title.String = 'Climate';
gridm
fname = strcat(figuredir,"climate.jpg");
print(fname,"-djpeg")


% 3. Print initial land cover maps
% --------------------------------
% (1) current (cleaned) land cover
figure(101);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
% h.Position = [1,31,1920,1093];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
cl = double(currentcleanlandcover);
cl(cl == 14) = 13;
cl(cl == 255) = NaN;
pcolorm(lats,lons,cl)
plotm(coastlat,coastlon,'Color','Black')
cclccm = [lightgreyval;reordcbp(1:12,:);reordcbp(14,:)];
colormap(cclccm)
c = colorbar;
c.Ticks = 0+13/28:13/14:14;
c.TickLabels = cellstr(["SNO/BAR/URB",biomes_igbp(2:13),biomes_igbp(15)]);
c.FontSize = 18;
ax = gca;
% ax.Position = [0.0156,0.11,0.8617,0.815];
ax.Position = [0,0.019,0.877,0.906];
ax.Title.FontSize = 24;
ax.Title.String = 'Current land cover';
gridm
fname = strcat(figuredir,"CurrentCleanLandCover.jpg");
print(fname,"-djpeg")

% (2) Current forests
figure(102);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
cf = double(currentforest);
cf(cf == 8) = 6;
cf(cf == 9) = 7;
cf(cf == 255) = NaN;
pcolorm(lats,lons,cf)
plotm(coastlat,coastlon,'Color','Black')
cfcm = [lightgreyval;reordcbp(1:5,:);reordcbp(8:9,:)];
colormap(cfcm)
c = colorbar;
c.Ticks = 0+7/16:7/8:7;
c.TickLabels = cellstr(["No forest",biomes_igbp(2:6),biomes_igbp(9:10)]);
c.FontSize = 18;
ax = gca;
ax.Position = [0,0.019,0.877,0.906];
ax.Title.FontSize = 24;
ax.Title.String = 'Current forest cover';
gridm
fname = strcat(figuredir,"CurrentForestCover.jpg");
print(fname,"-djpeg")

% (3) Current open lands
figure(103);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
ol = double(currentopenland);
ol(ol==7) = 1;
ol(ol==10) = 2;
ol(ol==12) = 3;
ol(ol==14) = 4;
ol(ol == 255) = NaN;
pcolorm(lats,lons,ol)
plotm(coastlat,coastlon,'Color','Black')
olcm = [lightgreyval;reordcbp(7,:);reordcbp(10,:);reordcbp(12,:);reordcbp(14,:)];
colormap(olcm)
c = colorbar;
c.Ticks = 0+4/10:4/5:4;
c.TickLabels = cellstr(["No open land",biomes_igbp([8,11,13,15])]);
c.FontSize = 18;
ax = gca;
ax.Title.FontSize = 24;
ax.Position = [0,0.019,0.877,0.906];
ax.Title.String = 'Current open lands';
gridm
fname = strcat(figuredir,"OpenLands.jpg");
print(fname,"-djpeg")


% 3. Print final land cover maps
% ------------------------------
% (1) Most likely forest
mlf = double(MLF);
mlf(mlf==u8noval) = NaN;
figure(104);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
pcolorm(lats,lons,mlf)
plotm(coastlat,coastlon,'Color','Black')
mlfcm = [lightgreyval;reordcbp(1:5,:)];
colormap(mlfcm)
c = colorbar;
c.Ticks = 0+5/12:5/6:5;
c.TickLabels = cellstr(["Not defined",biomes_igbp(2:6)]);
c.FontSize = 18;
ax = gca;
ax.Position = [0,0.019,0.877,0.906];
ax.Title.FontSize = 24;
ax.Title.String = 'Most likely forest';
gridm
fname = strcat(figuredir,"MostLikelyForest.jpg");
print(fname,"-djpeg")

% (2) Most likely open lands
mlol = double(MLOL);
mlol(mlol==u8noval) = NaN;
mlol(mlol==7) = 1;
mlol(mlol==10) = 2;
mlol(mlol==12) = 3;
mlol(mlol==14) = 4;
figure(110);clf
h = gcf;
h.Units = 'pixels';
h.Position = [1,208,1920,855];
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
load coastlines
pcolorm(lats,lons,mlol)
plotm(coastlat,coastlon,'Color','Black')
olcm = [lightgreyval;reordcbp(7,:);reordcbp(10,:);reordcbp(12,:);reordcbp(14,:)];
colormap(olcm)
c = colorbar;
c.Ticks = 0+4/10:4/5:4;
c.TickLabels = cellstr(["No open land",biomes_igbp([8,11,13,15])]);
c.FontSize = 18;
ax = gca;
ax.Position = [0,0.019,0.877,0.906];
ax.Title.FontSize = 24;
ax.Title.String = 'Most likely open lands';
gridm
fname = strcat(figuredir,"MostlikelyopenLands.jpg");
print(fname,"-djpeg")

% (3) Gridded forests
for gg = 1 : 8
    if gg == 1, var = PF01; vn = "PF01"; res = "0.1"; bx = "a)"; end
    if gg == 2, var = PF025; vn = "PF025"; res = "0.25"; bx = "b)"; end
    if gg == 3, var = PF05; vn = "PF05"; res = "0.5"; bx = "c)"; end
    if gg == 4, var = PF1; vn = "PF1"; res = "1"; bx = "d)"; end
    if gg == 5, var = PF25; vn = "PF25"; res = "2.5"; bx = "e)"; end
    if gg == 6, var = PF5; vn = "PF5"; res = "5"; bx = "f)"; end
    if gg == 7, var = PF10; vn = "PF10"; res = "10"; bx = "g)"; end
    if gg == 8, var = PF20; vn = "PF20"; res = "20"; bx = "h)"; end
    var = double(var);
    var(var==u8noval) = NaN;
    figure(10+gg);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,208,1920,865];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    pcolorm(lats,lons,var)
    plotm(coastlat,coastlon,'Color','Black')
    mlfcm = [lightgreyval;reordcbp(1:5,:)];
    colormap(mlfcm)
    c = colorbar;
    c.Ticks = 0+5/12:5/6:6;
    c.TickLabels = cellstr(["Not defined",biomes_igbp(2:6)]);
    c.FontSize = 36;
    ax = gca;
    ax.Position = [0.0156,0.0156,0.8146,0.899];
    ax.Title.FontSize = 48;
    ax.Title.String = strcat("Forest - ",res,"\circ x ",res,"\circ grid");
    annotation('textbox',[.03 .8 .1 .1],'String',bx,'FontSize',36,'EdgeColor',[1 1 1],...
        'FitBoxToText','on')
    gridm
    fname = strcat(figuredir,"GriddedForest_",vn,".jpg");
    print(fname,"-djpeg")
end

% (4) Gridded open lands
for gg = 1 : 8
    if gg == 1, var = POL01; vn = "POL01"; res = "0.1"; bx = "a)"; end
    if gg == 2, var = POL025; vn = "POL025"; res = "0.25"; bx = "b)"; end
    if gg == 3, var = POL05; vn = "POL05"; res = "0.5"; bx = "c)"; end
    if gg == 4, var = POL1; vn = "POL1"; res = "1"; bx = "d)"; end
    if gg == 5, var = POL25; vn = "POL25"; res = "2.5"; bx = "e)"; end
    if gg == 6, var = POL5; vn = "POL5"; res = "5"; bx = "f)"; end
    if gg == 7, var = POL10; vn = "POL10"; res = "10"; bx = "g)"; end
    if gg == 8, var = POL20; vn = "POL20"; res = "20"; bx = "h)"; end
    var = double(var);
    var(var==u8noval) = NaN;
    var(var==7) = 1;
    var(var==10) = 2;
    var(var==12) = 3;
    var(var==14) = 4;
    figure(20+gg);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [1,208,1920,865];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    pcolorm(lats,lons,var)
    plotm(coastlat,coastlon,'Color','Black')
    olcm = [lightgreyval;reordcbp(7,:);reordcbp(10,:);reordcbp(12,:);reordcbp(14,:)];
    colormap(olcm)
    c = colorbar;
    c.Ticks = 0+4/10:4/5:4;
    c.TickLabels = cellstr(["Not defined",biomes_igbp([8,11,13,15])]);
    c.FontSize = 36;
    ax = gca;
    ax.Position = [0.0156,0.0156,0.8146,0.899];
    ax.Title.FontSize = 48;
    ax.Title.String = strcat("Open Lands - ",res,"\circ x ",res,"\circ grid");
    annotation('textbox',[.03 .8 .1 .1],'String',bx,'FontSize',36,'EdgeColor',[1 1 1],...
        'FitBoxToText','on')
    gridm
    fname = strcat(figuredir,"GriddedOpenLands_",vn,".jpg");
    print(fname,"-djpeg")
end


% 3. Print result maps
% --------------------
valcat = [-1200,-800:200:-200,-100,-50,-20,20,50,100,200:200:800,1200];
nbcc = floor((length(valcat)-1)/2);
blueint = floor(lastblue/nbcc);
bluescale = jetmap(1:blueint:lastblue,:);
nlb = numel(bluescale)/3;
if nlb > nbcc, bluescale = bluescale([1:nlb-2,nlb],:); end
yellowint = floor((length(jetmap)-firstyellow)/nbcc);
yellowscale = jetmap(firstyellow:yellowint:length(jetmap),:);
nly = numel(yellowscale)/3;
if nly > nbcc, yellowscale = yellowscale([1,3:nly],:); end
catblind = cat(1,bluescale,greyval,yellowscale);
catblind = flip(catblind);

for ff = 1 : numel(outvarnames)
    
    eval(strcat("data = double(",outvarnames(ff),");"))
    
    data(data==i16noval) = NaN; %#ok<SAGROW>
    catdata = discretize(data,valcat);
    val1 = round(min(data,[],'all'));
    val2 = round(max(data,[],'all'));

    figure(ff);clf
    h = gcf;
    h.Units = 'pixels';
    h.Position = [8,157,1906,911];
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    pcolorm(lats,lons,catdata)
    plotm(coastlat,coastlon,'Color','Black')
    caxis([1 length(valcat)])
    colormap(catblind)
    c = colorbar;
    c.Ticks = 1:length(valcat);
    c.TickLabels = num2cell(valcat);
    c.FontSize = 20;
    c.Label.String = 'CO_2e [t CO_2 ha^-^1]';
    ax = gca;
    ax.Position = [0.0156,0.0156,0.8146,0.959];   
    ax.Title.FontSize = 24;
    ax.Title.String = maplegends(ff,1);
    bx = strcat("range is [",num2str(round(val1)),",",num2str(round(val2)),"]");
    leg = strcat("CO_2e of landcover change albedo-induced annual mean radiative forcing @ ToA",...
        " - kernels' ",maplegends(ff,2));
    annotation('textbox',[.702 .851 .121 .045],'String',bx,'FontSize',18,'EdgeColor','none',...
        'FitBoxToText','on')
    annotation('textbox',[.1578 .0073 .561 .055],'String',leg,'FontSize',18,'EdgeColor','none',...
        'FitBoxToText','on')
    annotation('textarrow',[.95 .95],[.2 .06],'String',"climate warming",'FontSize',18,...
        'TextRotation',90,'VerticalAlignment','top','HorizontalAlignment','center')
    annotation('textarrow',[.95 .95],[.78 .93],'String',"climate cooling",'FontSize',18,...
        'TextRotation',90,'VerticalAlignment','top','HorizontalAlignment','center')

%     
%     {maplegends(ff);...
%         strcat("\rm CO_2e of landcover change albedo-induced annual mean radiative forcing @ ToA ",...
%         " - range is [",num2str(round(val1)),",",num2str(round(val2)),"]")};
    gridm
    fname = strcat(figuredir,lower(outmapnames(ff)),".jpg");
    print(fname,"-djpeg")
    close(figure(ff))
end



% % (3) current land cover
% figure(100);clf
% h = gcf;
% h.Units = 'pixels';
% h.Position = [1,31,1920,1093];
% axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
% pcolorm(lats,lons,currentlandcover)
% plotm(coastlat,coastlon,'Color','Black')
% colormap(reordcbp)
% c = colorbar;
% c.Ticks = 1+15/32:15/16:16;
% c.TickLabels = cellstr(biomes_igbp(2:nbiomes));
% c.FontSize = 18;
% ax = gca;
% ax.Title.FontSize = 24;
% ax.Title.String = "Current land cover";
% gridm
% fname = strcat(figuredir,"CurrentLandCover.jpg");
% print(fname,"-djpeg")
% 


% % (6) Regions (when using only 9)
% reg = double(RE);
% reg(reg==u8noval) = NaN;
% figure(105);clf
% h = gcf;
% h.Units = 'pixels';
% h.Position = [1,31,1920,1093];
% axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
% load coastlines
% pcolorm(lats,lons,reg)
% plotm(coastlat,coastlon,'Color','Black')
% regcm = colorblindpalette([2:4,6,8:9,11,14:16],:);
% colormap(regcm)
% c = colorbar;
% c.Ticks = 1+10/22:9/10:10;
% c.TickLabels = cellstr(regionnames);
% c.FontSize = 18;
% ax = gca;
% ax.Title.FontSize = 24;
% ax.Title.String = 'Regions';
% gridm
% fname = strcat(figuredir,"Regions.jpg");
% print(fname,"-djpeg")

catblind = [0,0,0.5156;0,0,0.7188;0,0,0.9219;0,0.1250,1;0,0.3281,1;0,0.5312,1;0,0.9375,1;...
    0.7,0.7,0.7;1,1,0;1,0.5938,0;1,0.3906,0;1,0.1875,0;0.9844,0,0;0.7812,0,0;0.5781,0,0];
