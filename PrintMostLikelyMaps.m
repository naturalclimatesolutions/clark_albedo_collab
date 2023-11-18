% Print Most Likely Land Cover Maps (both Geotiffs and Matlab figures)

% Prepare parameters
% ******************
fileroot = strcat(regoutputfiles,"MostLikely005_");
indexfileroot = strcat(regoutputfiles,"ROinputs_");
inputdataformat = "uint8";
lowerresolution = 0.1;
method = "mode";
missinglandblock = 51;
mlindnames = ["openlandind","forestind","forwsaind","forwsaind"];

latlr = latlim005(2) - lowerresolution/2 : -lowerresolution : latlim005(1) + lowerresolution/2;
lonlr = lonlim005(1) + lowerresolution/2 : lowerresolution : lonlim005(2) - lowerresolution/2;
[lons,lats] = meshgrid(lonlr,latlr);


% Print maps
% **********
for ll = 1 : numel(mlnames)
    input = mlnames(ll);
    [lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,...
        indexfileroot);

    % Export Geotiff (high-resolution) map
    name = erase(mllongnames(ll)," ");
    fname = strcat(resfolder,name,"005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    lrvarname = strcat(mlnames(ll),"lowres");
    eval(strcat(lrvarname," = lowresmap;"));

    save(ReforOppfname,lrvarname{:},'-append')
    
    % Pring jpeg figure
    figure(ll); clf
    figuretype = "two panels";
    mapcolors = reordcbp;
    fname = strcat(figuredir,erase(mllongnames(ll)," "));
    albedofigure(lowresmap,IGBPAbbBio,lats,lons,fname,figuretype,mapcolors,...
        [],coastlat,coastlon," ",false);

    clear lowresmap highresmap
end


% Supplemental histogram
% **********************
biolabels = ["TrMF","TrDF","TrCF","TeMF","TeCF","BF","TrS","TeS","FlG","MGr","Tun","Med","Des","Man"];
semigeographicorder = [11,6,8,12,5,4,7,3,2,1,14,9,10,13];
o = semigeographicorder;
tables = ["MLOLdt","MLFCdt"];


ylim = round(ceil(max(sum(MLOLdt),sum(MLFCdt))*1.1),-1);
fignum = ["a.","b."];
for aa = 1 : numel(reforestationopp)
    atbl = areatable(:,:,aa);
    ctbl = TotalCO2byBiome_AObase(:,:,aa,ncont+1);
    stbl = JustCarbonbyBiome_AObase(:,:,aa,ncont+1);
    
    totco2 = sum(ctbl(o,:),2);
    posco2 = offsetcat(1:nocat)<50;
    totseq = sum(stbl(o,:),2);
    posco2val = sum(ctbl(o,posco2),2);
    
    figure(aa); clf
    h = gcf;
    h.Position = [500, 200, 977, 892];
    b = bar(atbl(o,:),'stacked');
    for ii = 1 : nocat
        %                 b(ii).FaceColor = 'flat';
        b(ii).FaceColor = aocolorbar(ii,:);
    end
    ax = gca;
    ax.FontSize = 20;
    if aa == 1
        ax.XTick = 1:nbiomes;
        ax.XTickLabel = BiomeAbb(o);
    else
        ax.XTick = [];
    end
    ax.YLim = [0 ylim];
    %             ax.TickDir = 'none';
    %             ax.YGrid = 'on';
    colormap(aocolorbar)
    c = colorbar;
    c.Ticks = 0:1/nocat:1;
    c.TickLabels = biolabels;
    c.Label.String = "Albedo Offset [%]";
    ax.YLabel.String = "Opportunity Area [Mha]";
    pos = ax.Position;
    xint = pos(3)/(nbiomes+1.4);
    xint2 = pos(3)/ (nbiomes+1.21);
    ysc = ax.YLim(2);
    bhgt = sum(atbl(o,:),2);
    b2hgt = bhgt - sum(atbl(o,posco2),2);
    b3hgt = bhgt - sum(atbl(o,posco2),2)/2;
    for bb = 1 : nbiomes
        bgx = pos(1) + (bb-0.3)*xint;
        bgx2 = pos(1) + (bb+0.5)*xint;
        bgx3 = pos(1) + .01 + (bb+0.32)*xint2;
        bgy = pos(2) + pos(4)/ysc*bhgt(bb);
        bgy2 = pos(2) + pos(4)/ysc*b2hgt(bb);
        bgy3 = pos(2) + pos(4)/ysc*b3hgt(bb);
        txt2 = strcat("(",num2str(round(totseq(bb),1)),")");
        annotation('textbox',[bgx bgy xint pos(4)/12],'String',...
            num2str(round(totco2(bb),1),'%.1f'),'EdgeColor','none','FontSize',18,...
            'FontWeight','bold','HorizontalAlignment','center')
        annotation('textbox',[bgx bgy xint pos(4)/20],'String',num2str(round(totseq(bb),1),...
            '%.1f'),'EdgeColor','none','FontSize',16,'HorizontalAlignment','center')
    end
    bgx = pos(1) + 0.2*xint;
    bgy = pos(2) +  0.9 * pos(4);
    annotation('textbox',[bgx bgy xint pos(4)/12],'String',fignum(aa),...
        'EdgeColor','none','FontSize',20,'FontWeight','bold','HorizontalAlignment','center')
    
    fname = strcat(figuredir,"Fig3",fignum(aa),reforestationopp(aa),"_AOarea_byBiome.jpg");
    print(fname,"-djpeg")
end
        
%     legend(flip(b)) % when using legend instead of colorbar but with the bins in the same order
%                       than on the graph







%% Additional maps for analysis
% ****************************

% for analysis, print modis
for yy = 1 : numel(modisyears)
    input = modisvars(yy);
    [~,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,cllist,mislist);
    
    eval(strcat(modisvars(yy),"lowres = lowresmap;"));

    fname = strcat(resfolder,input,"_005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);
end

% for analysis, print grids of MLF
for rr = 1 : numel(gridres)
    input = gridlist(rr,3);
    [~,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,...
        cllist,mislist,indexfileroot);

    fname = strcat("J:\ResultMaps\MostLikely\",input,"_005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    strcat("done with grid ",num2str(gridres(rr)))
end

% for analysis, print other interemediary maps
ll = 4;
intnames = ["mask";"ecoregion level (steps 1-3)";"local biome level (steps 1-4)";...
    "Added values from regional Biome/Climate (step 4)";...
    "Added values from global Biome/Climate (step 5)";'Added values from "other rules" (step 6)'];
for ii = 2 : numel(analysismaps)
    input = strcat(mlnames(ll),analysismaps(ii));
    [lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,...
        cllist,mislist,indexfileroot);

    name = strcat(erase(mllongnames(ll)," "),"_",analysismaps(ii));
    fname = strcat(resfolder,name,"005.tif");
    geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);

    
    eval(strcat("landcoverindex = ",mlindnames(ll),";"))
    nlc = numel(landcoverindex); extent = nlc+1;
    
    figure(ii); clf
    mlval = single(lowresmap);
    mlval(mlval==u8noval) = nan;
    for lc = 1 : nlc
        mlval(mlval==landcoverindex(lc)) = lc;
    end
    h = gcf;
    h.Units = 'pixels';
    h.Position = supercompmonitor;
    axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
    pcolorm(lats,lons,mlval)
    plotm(coastlat,coastlon,'Color','Black')
    mlfcm = [lightgreyval;reordcbp(landcoverindex,:)];
    colormap(mlfcm)
    c = colorbar;
    c.Ticks = 0+(extent-1)/(2*extent):(extent-1)/extent:nlc;
    c.TickLabels = cellstr(["Not defined";IGBPAbbBio(landcoverindex)]);
    c.FontSize = 18;
    ax = gca;
    ax.Position = scmpos;
    ax.Title.FontSize = 24;
    ax.Title.String = strcat(mllongnames(ll)," - ",intnames(ii));
    gridm
    fname = strcat(figuredir,name,".jpg");
    print(fname,"-djpeg")
    
    lrvarname = strcat(mlnames(ll),analysismaps(ii),"lowres");
    hrvarname = strcat(mlnames(ll),analysismaps(ii),"highres");
    eval(strcat(lrvarname," = lowresmap;"));
    eval(strcat(hrvarname," = highresmap;"));
    
    if exist(PotVegfname,'file') == 0
        save(PotVegfname,lrvarname{:},hrvarname{:})
    else
        save(PotVegfname,lrvarname{:},hrvarname{:},'-append')
    end
    clear lowresmap highresmap
end

% for analysis, print ecoregion rule
input = "forestsavannamap";
[lowresmap,highresmap] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
    rastersize005,blocksize005,lowerresolution,method,missinglandblock,...
    cllist,mislist,indexfileroot);

fname = strcat(resfolder,input,"005.tif");
geotiffwrite(fname,highresmap,R005,'TiffTags',cmptag);
extent = 2;

map = single(lowresmap);
map(map==255) = nan;
figure(100); clf
h = gcf;
h.Units = 'pixels';
h.Position = supercompmonitor;
axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
pcolorm(lats,lons,map)
plotm(coastlat,coastlon,'Color','Black')
mlfcm = reordcbp([3,8],:);
colormap(mlfcm)
c = colorbar;
c.Ticks = 1+(extent-1)/(2*extent):(extent-1)/extent:2;
c.TickLabels = cellstr(["Restricted to forests";"Woddy savannas allowed"]);
c.FontSize = 18;
ax = gca;
ax.Position = scmpos;
ax.Title.FontSize = 24;
ax.Title.String = strcat("Ecoregions potential forests' restrictions");
gridm
fname = strcat(figuredir,input,".jpg");
print(fname,"-djpeg")
    


