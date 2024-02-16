% NCS - Global Albedo - TNC-Clark Collaboration - Set Parameters

% Global maps and sub-regional dimensions
% ***************************************
% Maps dimensions
nlat1dg = 180;                          % 1 degree latitudes
nlon1dg = 360;                          % 1 degree longitudes
% 0.05 MODIS GCM scale, split into 20x20 degrees blocks
latlonscale = 0.05;                     % MODIS GCM scale
latlon05 = latlonscale;                 % Should call this when cleaning up ...
nlat05 = nlat1dg ./ latlon05;          % number of pixels in MODIS-GCM latitude-wise
nlon05 = nlon1dg ./ latlon05;          % number of pixels in MODIS-GCM longitude-wise
lat05 = 90 - latlon05/2 : -latlon05 : -90 + latlon05/2;
lon05 = -180 + latlon05/2 : latlon05 : 180 - latlon05/2;
biggestblock = 20;
blocksize05 = biggestblock / latlon05;
nblx = nlon05 / blocksize05;
nbly = nlat05 / blocksize05;
nblocks = nblx * nbly;
ncs05 = 90:-biggestblock:-90;
wcs05 = -180:biggestblock:180;
presland = true(nblocks,1);
ref05info = georasterinfo(reg05fname);
if ref05info.RasterReference.CellExtentInLatitude ~= latlon05, error("Wrong resolution"); end
rastersize05 = ref05info.RasterSize;
latlim05 = ref05info.RasterReference.LatitudeLimits;
lonlim05 = ref05info.RasterReference.LongitudeLimits;
firsti05 = NaN(nbly,1);
firstj05 = NaN(nblx,1);
sty = find(ncs05<=latlim05(2),1);
stx = find(wcs05<=lonlim05(1),1,'last');
edy = find(ncs05>latlim05(1),1,'last');
edx = find(wcs05<lonlim05(2),1,'last');
sti = 1; stj = 1;
for yy = sty : edy
    firsti05(yy) = sti;
    sti = sti + blocksize05;
end
for xx = stx : edx
    firstj05(xx) = stj;
    stj = stj + blocksize05;
end
clear sti stj yy xx sty stx edy edx
% 0.005 "Walker" scale, split into 10X10 degrees blocks (from 90N to -60S only)
bigblocks = 10;
latlon005 = 0.005;
blocksize005 = bigblocks / latlon005;
ncs = 90:-bigblocks:-90; nnbly = numel(ncs) - 1;
wcs = -180:bigblocks:180; nnblx = numel(wcs) - 1;
nbblocks = nnblx * nnbly;
block005description = zeros(nbblocks,5);
presbio = true(nbblocks,1);
preswlk = false(nbblocks,1);
firstiW005 = NaN(nnbly,1);
firstjW005 = NaN(nnblx,1);
walkerinfo = georasterinfo(WalkerBPotAGBfname);
if walkerinfo.RasterReference.CellExtentInLatitude ~= latlon005, error("Wrong resolution"); end
rastersize005 = walkerinfo.RasterSize;
latlim005 = walkerinfo.RasterReference.LatitudeLimits;
lonlim005 = walkerinfo.RasterReference.LongitudeLimits;
sty = find(ncs<=ceil(latlim005(2)),1);
stx = find(wcs<=lonlim005(1),1,'last');
edy = find(ncs>round(latlim005(1),4),1,'last');
edx = find(wcs<min(ceil(lonlim005(2)),180),1,'last');
sti = 1;stj = 1;
for yy = sty : edy
    firstiW005(yy) = sti;
    sti = sti + blocksize005;
end
for xx = stx : edx
    firstjW005(xx) = stj;
    stj = stj + blocksize005;
end
for yy = sty : edy
    for xx = stx : edx
        bb = xx + (yy-1) * nnblx;
        preswlk(bb) = true;
    end
end
for bb = 1 : nbblocks
    y = ceil(bb/nnblx);
    x = bb - (y-1)*nnblx;
    northlat = ncs(y); westlon = wcs(x);
    block005description(bb,:) = [bb,y,x,northlat,westlon];
end
clear sti stj yy xx sty stx edy edx bb x y northlat westlon
smtilesize = blocksize005 / 5;            % to use with parfor when needed
nsmtilex = blocksize005/smtilesize;
nsmtiley = blocksize005/smtilesize;
nsmtiles = nsmtilex * nsmtiley;
smtilesi = zeros(nsmtiles,2); smtilesj = zeros(nsmtiles,2);
for tl = 1 : nsmtiles
    ty = ceil(tl/nsmtilex);
    tx = tl - (ty-1)*nsmtilex;
    smtilesi(tl,:) = [(ty-1) * smtilesize + 1, ty * smtilesize];
    smtilesj(tl,:) = [(tx-1) * smtilesize + 1, tx * smtilesize];
end
clear tx ty tl

% Modis bounding coordinates (SIN grid)
opts = detectImportOptions(modisxlsfile,'Sheet','BoundingCoordinates');
modisccvars = opts.VariableNames;
modiscoord = readtable(modisxlsfile,opts,'Sheet','BoundingCoordinates'); clear opts



% Legends and variable names
% **************************
monthln = ["January";"February";...     % Months names
    "March";"April";"May";"June";"July";"August";"September";"October";"November";"December"];
dum = char(monthln);
monthnames = string(dum(:,1:3));        % months 3-letter names (used in some file names)
nmonths = length(monthln);              % number of months

albedo_types = ...                      % Albedo values found in albedo atlas
    ["snow-free black-sky","snow-free white-sky","snow-covered black-sky","snow-covered white-sky"];
albsn = ["nosnowblacksky",...           % Albedo short names (used as variable names)
    "nosnowwhitesky","snowcovblacksky","snowcovwhitesky"];
albid = cat(2,...                       % String values used in albedo file naming
    ["bsa";"wsa";"bsa";"wsa"],["";"";"snow_";"snow_"]);
nalb = length(albedo_types);            % number of albedo types

% 25 regions and continents
regionnames = ["Antarctica","Asiatic Russia","Australia/New Zealand","Caribbean","Central America",...
    "Central Asia","Eastern Africa","Eastern Asia","Eastern Europe","European Russia",...
    "Melanesia","Micronesia","Middle Africa","Northern Africa","Northern America",...
    "Northern Europe","Polynesia","South America","Southeastern Asia","Southern Africa",...
    "Southern Asia","Southern Europe","Western Africa","Western Asia","Western Europe"];
continents = ["Antarctica","Africa","Asia","Europe","North America","South America","Oceania"];
contregid = zeros(numel(regionnames),1);
for cc = 1 : numel(continents)
    cont = continents(cc);
    ci = contains(regionnames,cont);
    contregid(ci) = cc;
end
clear cc ci cont
contregid(ismember(regionnames,["Australia/New Zealand","Melanesia","Micronesia","Polynesia"])) = ...
    find(strcmp(continents,"Oceania"));
contregid(ismember(regionnames,["Caribbean","Central America","Northern America"])) = ...
    find(strcmp(continents,"North America"));

% IGBP - Land Cover
opts = detectImportOptions(IGBPLegendfile);
igbpLCvars = opts.VariableNames;
IGBPLegend = readtable(IGBPLegendfile,opts);
IGBPBiomes = string(table2array(IGBPLegend(:,strcmp(igbpLCvars,"LCnames"))));
IGBPAbbBio = string(table2array(IGBPLegend(:,strcmp(igbpLCvars,"LCAbb"))));
LC05ind = table2array(IGBPLegend(:,strcmp(igbpLCvars,"MCD12C1")));
LC005ind = table2array(IGBPLegend(:,strcmp(igbpLCvars,"MCD12Q1")));
forwsabiomes = IGBPBiomes(contains(IGBPBiomes,["forest","Woody"]));
forestsavannabiomes = IGBPAbbBio(contains(IGBPBiomes,["forest","savanna"],"IgnoreCase",true));
forwsaind = LC005ind(contains(IGBPBiomes,["forest","Woody"]));
forsavind = LC005ind(contains(IGBPBiomes,["forest","savanna"],"IgnoreCase",true));
allsavind = LC005ind(contains(IGBPBiomes,"savanna","IgnoreCase",true));
wsaind = LC005ind(contains(IGBPBiomes,"Woody"));
forestbiomes = IGBPBiomes(contains(IGBPBiomes,"forest"));
forestind = LC005ind(contains(IGBPBiomes,"forest"));
openlandbiomes = IGBPBiomes(contains(IGBPBiomes,["Open","Cropland","Grassland"]));
openlandind = LC005ind(contains(IGBPBiomes,["Open","Cropland","Grassland"]));
landind = LC005ind(~contains(IGBPBiomes,"Water"));
waterind = LC005ind(contains(IGBPBiomes,"Water"));

% Dinerstein & al (2017) Ecoregions
opts = detectImportOptions(ecoregionxlsfile,'Sheet','Ecoregions2017');
ecovarnames = opts.VariableNames;
ecoregions = readtable(ecoregionxlsfile,opts,'Sheet','Ecoregions2017');
ecoindex = table2array(ecoregions(:,strcmp(ecovarnames,"OBJECTID")));
ecoid = table2array(ecoregions(:,strcmp(ecovarnames,"ECO_ID")));
realeconames = string(table2array(ecoregions(:,strcmp(ecovarnames,"ECO_NAME"))));
econames = erase(realeconames," ");
biomeid = table2array(ecoregions(:,strcmp(ecovarnames,"BIOME_NUM")));
biomelist = unique(biomeid); nbiomes = numel(biomelist);
dum = string(table2array(ecoregions(:,strcmp(ecovarnames,"BIOME_NAME"))));
biomenames = strings(nbiomes,1);
ecoperbiome = zeros(nbiomes,numel(ecoid)); mec = 0; necoperbiome = zeros(nbiomes,1);
for bb = 1 : nbiomes
    biomenames(bb) = dum(find(biomeid==biomelist(bb),1));
    ecobb = sort(ecoid(biomeid==bb));
    neco = numel(ecobb);
    necoperbiome(bb) = neco;
    ecoperbiome(bb,1:neco) = ecobb;
    if numel(ecobb) > mec, mec = numel(ecobb); end
end
tundradesertsindex = biomelist(contains(biomenames,["Tundra","Deserts"]));
forestedecoregions = ecoid(contains(realeconames,["forest","varzea","taiga","mangrove","pantanos",...
    "yungas","piney","pine barrens"],'IgnoreCase',true));
ecoperbiome = ecoperbiome(:,1:mec);
clear opts dum bb mec ecobb neco

% Koppen-Geiger climate zones
opts = detectImportOptions(KGLegend);
kglvarnames = opts.VariableNames;
kglegend = readtable(KGLegend,opts);
climateid = table2array(kglegend(:,strcmp(kglvarnames,"ID")));
kglegnames = table2array(kglegend(:,strcmp(kglvarnames,"Legend")));
polardesertindex = climateid(contains(kglegnames,["desert","Polar"]));
clear opts

% Walker classes
WlkClassID = 0:7;
WlkClassDef = ["NoData","Nonwoody", "R/L", "MM/L", "M/L", "R/H", "MM/H", "M/H"];
Wsc = WlkClassID(contains(WlkClassDef,"R/"));

% Kernels' specific names
kernels = ["CAM3","CAM5",...            % Radiative kernel names
    "ECHAM6","HADGEM2","CACK1","HADGEM3"];
nkernels = numel(kernels);              % number of radiative kernels
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
            kernelfilenames(kk) = "HADGEM3\HadGEM3-GA7.1_TOA_kernel_L85.nc";
            kernelvarnames(kk) = "albedo_sw";
            kernelscale(kk) = 100;
        otherwise
            error("This kernel has not been defined yet!")
    end
end
clear kk kname


% GeoTiff file export parameters
% ******************************
cmptag.Compression = 'LZW';             % LZW compression is better than the default ...
% a. 0.005 scale and [-60 90] latitudes
R005 = georefcells(latlim005,lonlim005,rastersize005,'ColumnsStartFrom','north');
% b. 0.05 scale and [-90 90] latitudes
R05 = georefcells(latlim05,lonlim05,rastersize05,'ColumnsStartFrom','north');


% File lists
% **********
dum = deblank(string(ls(biomassfolder)));
biofilelist = dum(endsWith(dum,"tif")); clear dum
        
dum = deblank(string(ls(modis0008folder)));
modisfilelist0008 = dum(endsWith(dum,"tif")); clear dum
dum = deblank(string(ls(modisorigdata)));
modisoriginalfilelist = dum(endsWith(dum,"hdf")); clear dum
dum = deblank(string(ls(strcat(modis005folder,"2001\"))));
modis2001filelist005 = dum(endsWith(dum,"tif")); clear dum
dum = deblank(string(ls(strcat(modis005folder,"2010\"))));
modis2010filelist005 = dum(endsWith(dum,"tif")); clear dum

dum = deblank(string(ls(ecofolder)));
ecofilelist = dum(endsWith(dum,"tif")); clear dum


% Other variables and parameters
% ******************************
cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double","logical"];
mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan,false];
u16noval = mislist(cllist=="uint16");
u8noval = mislist(cllist=="uint8");
i16noval = mislist(cllist=="int16");

despct = [50,60,75:5:95];
mappct = [5,10,25,50,75,90,95];
ostats = ["mean","std","median","25%prctile","75%prctile"];
gridres = [0.01,0.025,0.05,0.1,0.25,0.5,1,2.5,5,10];
kerstat = ["median","min","max"];

Ca_PgC  = 400*2.13;                     % current base CO2 in atmos [Pg C]
Aglobe  = 5.1007e14;                    % global surface area [m2]
nodata_albedo = 2^15-1;                 % Albedo "nodata" value
scale_albedo = 0.001;                   % Albedo scale value
dm2C = 0.47;                            % Dry biomass to Carbon conversion (IPCC)
GWP100 = 1.9;                           % GWP*100 for albedo to account for land/ocean response
earthellipsoid = referenceEllipsoid('wgs84','m');
albedostats = ["median","min","max"];
reforestationopp = ["Walker","Griscom","Bastin"];
allreforestationopp = cat(2,reforestationopp,"CombinedOpp","Globally");
offsetcat = ...                         % Albedo Offset categories 
    [-10000,25,50,75,100,10000];  %  (we now reversed the sign to avoid double-negation)
aothreshold = 50;
seqcat = [-2000,-800:200:-200,-100,-40,40,100,200:200:800,2000];
% seqcat = [-2000,-400,-200,-100,-40,40,100,200,400,2000];
nocat = numel(offsetcat) - 1;
nscat = numel(seqcat) - 1;
landcoverpriority = "newest";
prevalentlandcoverthreshold = 1;        % precent area below which land cover is ignored in
                                        % "most likely" calculation
minecopixels = 20;                      % Same idea, but as a minimum number of pixels (as the 1% 
minbiomepixels = 100;                   % was too restrictive at the ecoregion/biome level)

singlepathstatistics = ...              % list of statistics to be performed over kernels
    ["min","max","mean","median","range","standard deviation"];
nstats = numel(singlepathstatistics);



% Figures parameters
% ******************
load coastlines;
colorblindpalette = [0,0,0;0,73,73;0,146,146;255,109,182;255,182,219;73,0,146;0,109,219;...
    182,109,255;109,182,255;182,219,255;146,0,0;146,73,0;219,109,0;36,255,36;255,255,219;...
    255,255,109]/255;
landcoverpalette = colorblindpalette([1:3,9:10,6:7,4:5,14,11,8,12:13,15:16],:);
reordcbp = colorblindpalette([1:3,9:10,6:7,4:5,14,11,8,12:13,15:16],:);
lightgreyval = [.9,.9,.9];
greyval = [0.7,0.7,0.7];
BCredtoblue11 = [165,0,38;215,48,39;244,109,67;253,174,97;254,224,144;255,255,191;...
    224,243,248;171,217,233;116,173,209;69,117,180;49,54,149]/255;
BCpurpletogreen8 = [118,42,131;153,112,171;194,165,207;231,212,232;217,240,211;166,219,160;90,174,97;27,120,55]/255;
BCpurpletogreen6 = [118,42,131;175,141,195;231,212,232;217,240,211;127,191,123;27,120,55]/255;
BCpurpletogreen4 = [123,50,148;194,165,207;166,219,160;0,136,55]/255;
BCbrowntoteal6 = [140,81,10;216,179,101;246,232,195;199,234,229;90,180,172;1,102,94]/255;
aocolorbar = flip(cat(1,BCpurpletogreen6(1:3,:),BCpurpletogreen4(3:4,:)));
aosenscolorbar = cat(1,lightgreyval,BCbrowntoteal6([1,2,3,5],:));
co2colorbar = cat(1,[0.5312 0 0],BCredtoblue11,[0,0,0.5156]);
co2colorbaror = redtobluecolorbar(seqcat,7,true,"redtoblue",false);
laptopscreen = [1 31 1536 758];
% homemonitor = [240 895 1680 944];
homemonitor = [16 93 1654 850];
workmonitor = [1 31 1920 1093];
supercompmonitor = [1 31 2560 1333];
newsupcmpmonitor = [256 271 2047 1026];
ltppos = [0.019 0.019 0.8 0.95];
lptmpos = [0.0426 0.019 0.8128 0.95];
scmpos = [0.0426, 0.0190, 0.8344, 0.9060];
scmposwm = [0.0426    0.0190    0.8194    0.9500];
scmantpos = [0.0426, 0.0190, 0.8344, 0.95];


% Save parameters
% ***************
save(parameterfile)
