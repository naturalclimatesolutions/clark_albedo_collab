% Export Sensitivity map
% **********************

indexfileroot = strcat(regoutputfiles,"ROinputs_");
fileroot = strcat(regoutputfiles,"ROstats005_");
inputdataformat = "uint8";
lowerresolution = 0.1;
missinglandblock = 51;

latlr = latlim005(2) - lowerresolution/2 : -lowerresolution : latlim005(1) + lowerresolution/2;
lonlr = lonlim005(1) + lowerresolution/2 : lowerresolution : lonlim005(2) - lowerresolution/2;
[lons,lats] = meshgrid(lonlr,latlr);

[AOsenslowres,AOsenshighres] = printmapprep('AOsensmap',fileroot,preswlk,latlon005,inputdataformat,...
    rastersize005,blocksize005,lowerresolution,'mode',missinglandblock,indexfileroot);

% Export Geotiff (high-resolution) map
fname = strcat(resfolder,"AOSensitivity005.tif");
geotiffwrite(fname,AOsenshighres,R005,'TiffTags',cmptag);

save(ReforOppfname,"AOsenslowres","AOsenshighres",'-append')
load(ReforOppfname,'sensarreas50')

% Print figure
pval = round(sensarreas50 ./ sum(sensarreas50) * 100);
labels = [strcat("Not Defined \fontsize{16} \it \newline (",num2str(pval(1)),"% land area)");...
    strcat("always >50%AO \fontsize{16} \it \newline (",num2str(pval(2)),"% land area)");...
    strcat(">50% in max&med \newline albedo effect \fontsize{16} \it \newline (",num2str(pval(3)),"% land area)");
    strcat(">50% only in max \newline albedo effect \fontsize{16} \it \newline (",num2str(pval(4)),"% land area)");
    strcat("always <50%AO \fontsize{16} \it \newline (",num2str(pval(5)),"% land area)")];
figure(100);clf
figurename = strcat(figuredir,"AOSensitivity");
thisfiguretype = "two panels";
axislegend = "Albedo Offset [%]";
h = albedofigure(AOsenslowres,sensitivitycategories,lats,lons,figurename,thisfiguretype,...
    aosenscolorbar,0,coastlat,coastlon,axislegend,false,labels);



% Just checking
% *************
fileroot = strcat(regoutputfiles,"RadForcing005_");
inputdataformat = "single";
method = "mean";
hiresnames = ["AlbedoRadiativeForcing","AlbedoRFmin","AlbedoRFmax"];
albedostatname = 

for ll = 1 : bdnames
    input = strcat("AO",bdnames(ll));
    [lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot);
    
    lowresmap(lowresmap<offsetcat(1)) = offsetcat(1);
    lowresmap(lowresmap>offsetcat(numel(offsetcat))) = offsetcat(numel(offsetcat));
    highresmap(highresmap<offsetcat(1)) = offsetcat(1);
    highresmap(highresmap>offsetcat(numel(offsetcat))) = offsetcat(numel(offsetcat));

    % Export Geotiff (high-resolution) map
    fname = strcat(resultmaps,hiresnames(ll),"_005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    lrvarname = strcat(input,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));

    save(ReforOppfname,lrvarname,'-append')
    
    if contains(input,"AO")
        thiscolorbar = aocolorbar;
        thesecategories = offsetcat;
        axislegend = "Albedo Offset [%]";
    else
        thiscolorbar = co2colorbar;
        thesecategories = seqcat;
        axislegend = 'Mg CO_2e ha^-^1';
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






