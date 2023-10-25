% Print Albedo Carbon Maps
% ************************

% Prepare parameters
% ******************
fileroot = strcat(regoutputfiles,"RadForcing005_");
indexfileroot = strcat(regoutputfiles,"ROinputs_");
inputdataformat = "single";
lowerresolution = 0.1;
method = "mean";
missinglandblock = 51;
maps = ["RFmed005","WlkTCO2RF","NCIbase","AObase"];
oppmaps = strcat(reforestationopp,"NCI");
figurenumber = ["Fig1a_","Fig1b_","Fig1c_","FigS1_"];
hiresnames = ["AlbedoRadiativeForcing","WalkerCO2RFlim","NetClimateImpact","AlbedoOffset"];
statname = ["AlbedoRF_AreaCountStats";"CarbonOnlyRFlimited_AreaCountStats";"GloballyNCI_AreaCountStats";...
    "AlbedoOffset_AreaCountStats"];
oppstatname = strcat(oppmaps,"_AreaCountStats");
hrvarname = strings(numel(maps),1);
oppfignumber = ["FigS3c_","FigS3a_","FigS3b_"];
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
        thiscolorbar = aocolorbar;
        thesecategories = offsetcat;
        axislegend = "Albedo Offset [%]";
        arrows = false;
    else
        thiscolorbar = co2colorbar;
        thesecategories = seqcat;
        axislegend = 'Mg CO_2e ha^-^1';
        arrows = true;
    end
    
    eval(strcat("prctvals =",statname(ll),";"))
    
    % Pring jpeg figure
    figurename = strcat(figuredir,figurenumber(ll),hiresnames(ll),".jpeg");
    figure(ll); clf
    h = albedofigure(lowresmap,thesecategories,lats,lons,figurename,positions,prctvals,...
    coastlat,coastlon,axislegend,thiscolorbar,arrows);
    
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
    end
    [lowresmap,~] = printmapprep("NCIbase",fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot,mask);


    lrvarname = strcat(oppname,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));
    save(ReforOppfname,lrvarname{:},'-append')
    
    thiscolorbar = co2colorbar;
    thesecategories = seqcat;
    axislegend = 'Mg CO_2_e ha^-^1';
    
    eval(strcat("prctvals =",oppstatname(ll),";"))
    
    % Pring jpeg figure
    figurename = strcat(figuredir,oppfignumber(ll),oppname,".jpeg");
    figure(ll); clf
    h = albedofigure(lowresmap,thesecategories,lats,lons,figurename,positions,prctvals,...
    coastlat,coastlon,axislegend,thiscolorbar,true);
    
    clear lowresmap highresmap
end