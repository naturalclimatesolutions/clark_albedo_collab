% Print Albedo Carbon Maps
% ************************
% I now print them minimally and import them in illustrator to finalize

% Prepare parameters
% ******************
fileroot = strcat(regoutputfiles,"RadForcing005_");
indexfileroot = strcat(regoutputfiles,"ROinputs_");
inputdataformat = "single";
lowerresolution = 0.1;
method = "mean";
missinglandblock = 51;
maps = ["GWPmed005","WlkTCO2RF","NCIbase","AObase"];
oppmaps = cat(2,strcat(reforestationopp,"NCI"),"CombinedOppNCI");
figurenumber = ["FigS1a_","FigS1b_","FigS1c_","FigS1d_"];
hiresnames = ["AlbedoRadiativeForcing","WalkerCO2RFlim","NetClimateImpact","AlbedoOffset"];
hiresoppmaps = ["AlbedoRadiativeForcing","WalkerCO2RFlim","NetClimateImpact","AlbedoOffset"];
statname = ["AlbedoGWP_AreaCountStats";"CarbonOnlyRFlimited_AreaCountStats";"GloballyNCI_AreaCountStats";...
    "AlbedoOffset_AreaCountStats"];
hrvarname = strings(numel(maps),1);
oppfignumber = ["Fig2e_","Fig2a_","Fig2c_",""];
axeslabels = 0;
arrows = false;
load(ReforOppfname, statsarraylist{:})


latlr = latlim005(2) - lowerresolution/2 : -lowerresolution : latlim005(1) + lowerresolution/2;
lonlr = lonlim005(1) + lowerresolution/2 : lowerresolution : lonlim005(2) - lowerresolution/2;
[lons,lats] = meshgrid(lonlr,latlr);
switch computer
    case '7264'
        positions.figure = workmonitor;
        positions.axes = scmposwm;
    case '6381'
        positions.figure = supercompmonitor;
        positions.axes = scmantpos;
    case '8904'
        positions.figure = [25 15 13.9960 4.6700];
%         positions.figure = [239 295 2047 683];
        positions.axes = [0.0426 0.019  0.8194 0.95];
end


% Print global maps
% *****************
for ll = 1 : numel(maps)
    input = maps(ll);
    [lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot);
    
    % Export Geotiff (high-resolution) map
    fname = strcat(resultmaps,hiresnames(ll),"_005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    lrvarname = strcat(input,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));

    save(ReforOppfname,lrvarname{:},'-append')
    
    if contains(input,"AO")
        lowresmap(lowresmap<offsetcat(1)) = offsetcat(1);
        lowresmap(lowresmap>offsetcat(numel(offsetcat))) = offsetcat(numel(offsetcat));
        thiscolorbar = aocolorbar;
        thesecategories = offsetcat;
        axislegend = "Albedo Offset [%]";
%         arrows = false;
    else
        thiscolorbar = co2colorbar;
        thesecategories = seqcat;
        axislegend = 'Mg CO_2e ha^-^1';
%         arrows = true;
    end

    thisfiguretype = "two panels";
    eval(strcat("prctvals =",statname(ll),";"))
    
    % Pring figure
    figurename = strcat(figuredir,figurenumber(ll),hiresnames(ll));
    figure(ll); clf
    h = albedofigure(lowresmap,thesecategories,lats,lons,figurename,thisfiguretype,thiscolorbar,...
    prctvals,coastlat,coastlon,axislegend,false);

    clear lowresmap highresmap
end

% Print only high res for data availability
% *****************************************
invars = ["RFmin005","RFmax005","NCIwlk85"];
highresnames = ["AlbedoRadiativeForcing_min.tif","AlbedoRadiativeForcing_max.tif",...
    "ESATruncNetClimateImpact.tif"];
for ll = 1 : numel(invars)
    input = invars(ll);
    [lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot);
    
    % Export Geotiff (high-resolution) map
    fname = strcat(resultmaps,highresnames(ll));
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    lrvarname = strcat(input,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));

    save(ReforOppfname,lrvarname{:},'-append')
end


% Print opportunity NCI maps
% **************************
for ll = 1 : numel(oppmaps)
    oppname = oppmaps(ll);
    switch oppname
        case "WalkerNCI"
            mask = "walkeroppcat";
        case "GriscomNCI"
            mask = "griscom";
        case "BastinNCI"
            mask = "selectedbastin";
        case "CombinedOppNCI"
            mask = "combinedopp";
    end
    [lowresmap,highresmap] = printmapprep("NCIbase",fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot,mask);


    % Export Geotiff (high-resolution) map
    fname = strcat(resultmaps,oppname,"_005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    lrvarname = strcat(oppname,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));
    save(ReforOppfname,lrvarname{:},'-append')
    
    thiscolorbar = co2colorbar;
    thesecategories = seqcat;
    thisfiguretype = "six panels";
    axislegend = 'Mg CO_2e ha^-^1';
    
    eval(strcat("prctvals =",oppstatname(ll),";"))
    
    % Pring jpeg figure
    figurename = strcat(figuredir,oppfignumber(ll),oppname);
    figure(ll); clf
    h = albedofigure(lowresmap,thesecategories,lats,lons,figurename,thisfiguretype,thiscolorbar,...
    prctvals,coastlat,coastlon,axislegend,false);
    
    clear lowresmap highresmap
end