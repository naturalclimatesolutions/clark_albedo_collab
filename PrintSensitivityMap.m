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

for kk = 1 : bdnames
    eval(strcat("input = AO",bdnames(kk)));
    




