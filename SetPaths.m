% NCS - Global Albedo - TNC-Clark Collaboration - Set File & Folder Paths

% Computer
% ********
computer = "8904";                      % either 7264 (laptop) or 6381 (supercomputer)
switch computer
    case "8916"
        dataroot = "F:\";
        homeroot = "C:\Work\";
        intmat = strcat(dataroot,"IntermediaryFiles\MATfiles\");
    case "6381"
        dataroot = "J:\";
        homeroot = "D:\";
        intmat = strcat(dataroot,"IntermediaryFiles\MATfiles\");
    case "8904"
        dataroot = "F:\";
        homeroot = "D:\";
        intmat = "E:\GlobalAlbedo\MATFiles\";

end

% Files and Folders
% *****************
% Folders
origdatafolder = strcat(dataroot,"OriginalDatasets\");
intdatafolder = strcat(dataroot,"IntermediaryFiles\");
inputdatafolder = strcat(dataroot,"IntermediaryFiles\GlobalAlbedoInputs\");
matfolder = strcat(homeroot,"GlobalAlbedo\MatlabFiles\");
resfolder = strcat(homeroot,"GlobalAlbedo\Results\");
boxpath = "C:\Users\nhasler\Box\Albedo (Clark-TNC)";

% Data in MODIS GCM 0.05 degrees
MODISSnowFilesFolderPath = "G:\OriginalDatasets\SnowData\MODIS_MOD10CM_v6.1\";
ResultFolderPath = "G:\ProcessedDatasets\SnowCover\";
landcovermap20 = strcat(origdatafolder,"MODISLandCover\MCD12C1.A2020001.006.2021362215328.hdf");
landcovermap01 = strcat(origdatafolder,"MODISLandCover\MCD12C1.A2001001.006.2018053185512.hdf");
albedosnowfree = strcat(origdatafolder,"AlbedoAtlas\snow_free_hierarchical_v3\");
albedosnowcovered = strcat(origdatafolder,"AlbedoAtlas\snow_covered_hierarchical_v3\");
albedoheader = strcat(albedosnowcovered,"hierarchical_snow_albedo_igbp_0.05.bsa_shortwave.Jan.hdr");
albedopaths = cat(1,albedosnowfree,albedosnowfree,albedosnowcovered,albedosnowcovered);
snowcoverdata = strcat(origdatafolder,"SnowData\MODIS_MOD10CM_v6.1\");
climsnowcover = strcat(intdatafolder,"SnowCover\");
solarradiation = strcat(origdatafolder,"SolarFlux\");
radiativekernel = strcat(origdatafolder,"RadiativeKernels\");
reg05fname = strcat(inputdatafolder,"Regions05.tif");
ecoregion05fname = strcat(inputdatafolder,"Ecoregions05.tif");

% Data in 0.005 degrees (WGS84)
ecofolder = strcat(intdatafolder,"Ecoregions10x10Tiles\");
modis005folder = strcat(intdatafolder,"MODISLandCoverWGS84_005\");
modisyears = ["2001","2010"];
walkerfolder = strcat(intdatafolder,"Walker\");
WalkerBPotAGBfname = strcat(walkerfolder,"Walker_BPotAGB_WGS84.tif");
WalkerBPotBGBfname = strcat(walkerfolder,"Walker_BPotBGB_WGS84.tif");
WalkerBPotLowPctfname = strcat(walkerfolder,"Walker_BPotLowValAGB_WGS84.tif");
WalkerBPotHighPctfname = strcat(walkerfolder,"Walker_BPotHighValAGB_WGS84.tif");
ecoregion005fname = strcat(inputdatafolder,"Ecoregions005.tif");
reg005fname = strcat(inputdatafolder,"WorldRegions005.tif");
KGClimatefname = strcat(inputdatafolder,"KoppenGeiger005.tif");
matecobiofolder = strcat(intdatafolder,"MATfiles\EcoregionBiomass\");
subdata10 = strcat(matecobiofolder,"SubRegionData\");
subdataeco = strcat(matecobiofolder,"EcoregionFiles\");
if ~exist(subdataeco,'dir'), mkdir(subdataeco); end
Bastinfname = strcat(inputdatafolder,"BastinRestorationPotential005.tif");
BastinTPfname = strcat(inputdatafolder,"BastinTotalPotential005.tif");
Griscomfname = strcat(inputdatafolder,"GriscomNCSWGS84005.tif");

% Data in 0.0008 resolution (WGS84)
modis0008folder = strcat(intdatafolder,"MODISLandCoverWGS84_00088\");
biomassfolder = strcat(origdatafolder,"ESACCIBiomass\2010\");

% Other maps
modisorigdata = strcat(origdatafolder,"MODISLandCover\MCD12Q1_2010\");

% Excel files
modisxlsfile = strcat(origdatafolder,"MODISTiles\MODISBoundingCoordinates.xlsx");
ecoregionxlsfile = strcat(inputdatafolder,"Ecoregions2017.xlsx");
KGLegend = strcat(origdatafolder,"Koppen-Geiger_Beck\legend.xlsx");
IGBPLegendfile = strcat(strcat(origdatafolder,"MODISLandCover\IGBPLegend.xlsx"));

% Output MATLAB files
parameterfile = strcat(matfolder,"GAParameters.mat");
sparameterfile = strcat(matfolder,"GASecondaryParameters.mat");
pathsfile = strcat(matfolder,"GAlocalpaths.mat");
regoutputfiles = strcat(intmat,"SubData\");
resultsGlasgowCOP = strcat(dataroot,"ResultMaps\GlasgowCOP\");
resultmaps = strcat(resfolder,"Maps\");
EgyptCOP = strcat(intmat,"EgyptCOP.mat");
AlbSubfolder = strcat(homeroot,"GlobalAlbedo\SubData\");
ForestSavannaFile = strcat(intmat,"ForestSavanna.mat");
ReforOppfnameold = strcat(intmat,"ReforestationOpportunities.mat");
EcoStatfname = strcat(matfolder,"EcoregionStatistics.mat");
ReforOppfname = strcat(matfolder,"ReforestationOpportunities.mat");
ReforOppLAfname = strcat(matfolder,"ReforestationOpportunitiesHighRes.mat");
PotVegfname = strcat(matfolder,"PotentialVegetation.mat");
W2Wfname = strcat(intmat,"WalltoWall.mat");
ESACCIfname = strcat(intmat,"ESACCI.mat");
figuredir = strcat(matfolder,"Figures\");
if ~exist(figuredir,'dir'), mkdir(figuredir); end
analysisfiguredir = strcat(intdatafolder,"Figures\");
if ~exist(analysisfiguredir,'dir'), mkdir(analysisfiguredir); end

save(pathsfile)