% NCS - Global Albedo - TNC-Clark Collaboration - Reforestation Opportunity Paper
%
% This code is a "master" code, where I try to condense and clean-up all the bits and pieces of code
%   produced along the way to complete the Reforestation Opportunities paper (along with the Glasgow
%   COP set of data, referenced in that paper)
%
% Created 12/16/2022 by Natalia Hasler, Clark University


clearvars

% Define paths/parameters
% ***********************
SetPaths
SetParameters

% Import base datasets (atlas, kernels, land-cover, snow conditions, solar radiation)
% ********************
generate_snow_cover_climatology
ImportBaseData

% Calculate Albedo change-induced climate forcing for specific land-covers (GlasgowCOP - data availability)
% ************************************************************************
ComputeSinglePathwayRF
ExportSinglePathwayMaps

% Import (or calculate) Maps @ 0.005 resolution
% *********************************************
ImportROdata
AdditionalROData

% Generate Most Likely Land Cover Maps
% ************************************
EcoregionStatistics
GenerateMostLikelyLandCover
PrintMostLikelyMaps

% Calculate Albedo-change induced climate forcing for most likely transitions
% ***************************************************************************
MLAlbedoInducedForcing
CombineAlbedoCarbon

% Calculate area statistics for base and opportunity maps
% *******************************************************
CalculateMapAreaStatistics
PrintAlbedoCarbonMaps

% Opportunity Stats and Albedo Sensitivity (now we limit it to albedo only)
% ****************************************
CalculateOpportunityStats
ExportOpportunityTables
PrintOpportunityBarCharts
PrintSensitivityMap

% Other carbon map (at this point, we only use the ESA-CCI-truncated Walker)
% ****************
ImportESACCIBiomass
CombineESAWalkerAlbedo
ExportESAWalkerAlbedo





