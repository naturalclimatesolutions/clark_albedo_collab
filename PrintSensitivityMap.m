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

% Print figure
pval = round(sensarreas50 ./ sum(sensarreas50) * 100);
labels = [strcat("Not Defined \fontsize{16} \it \newline (",num2str(pval(1)),"% land area)");...
    strcat("always <50%AO \fontsize{16} \it \newline (",num2str(pval(2)),"% land area)");...
    strcat(">50% in max&med \newline albedo effect \fontsize{16} \it \newline (",num2str(pval(3)),"% land area)");
    strcat(">50% only in max \newline albedo effect \fontsize{16} \it \newline (",num2str(pval(4)),"% land area)");
    strcat("always >50%AO \fontsize{16} \it \newline (",num2str(pval(5)),"% land area)")];
thesecategories = 0:4; nlc = numel(thesecategories); extent =nlc + 1;
axislegend = "Albedo Offset [%]";
figure(10); clf
mlval = single(AOsenslowres);
mlval(mlval==u8noval) = nan;
h = gcf;
h.PaperOrientation = 'landscape';
h.Units = 'pixels';
h.Position = laptopscreen;
    axesm('MapProjection','robinson','Frame','on','MLineLocation',20,'PLineLocation',20,...
        'Grid','on','MeridianLabel','on','MLabelLocation',60,'ParallelLabel','on',...
        'PLabelLocation',20,'LabelRotation','on','FontSize',18);
pcolorm(lats,lons,mlval)
plotm(coastlat,coastlon,'Color','Black')
colormap(aosenscolorbar)
c = colorbar;
c.Ticks = nlc/(2*extent):(nlc-1)/nlc:nlc;
c.TickLabels = cellstr(labels);
c.FontSize = 18;
ax = gca;
ax.Position = ltppos;
fname = strcat(figuredir,"AOSensitivity");
print(strcat(fname,".jpg"),"-djpeg")
print(strcat(fname,".pdf"),"-dpdf",'-bestfit')


