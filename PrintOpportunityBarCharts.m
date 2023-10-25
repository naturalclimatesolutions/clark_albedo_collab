% Plot bar charts (median values only)
% ***************
BiomeAbb = ["TrMF","TrDF","TrCF","TeMF","TeCF","BF","TrS","TeS","FlG","MGr","Tun","Med","Des","Man"];
semigeographicorder = [11,6,8,12,5,4,7,3,2,1,14,9,10,13];
aolabels = num2cell(offsetcat);
aolabels(1) = num2cell(strcat("< ",num2str(offsetcat(2))));
aolabels(numel(offsetcat)) = num2cell(strcat("> ",num2str(offsetcat(numel(offsetcat)-1))));

areatable = AreasbyBiome_AObase(:,:,1:numel(reforestationopp),ncont+1);
ylim = round(ceil(max(sum(areatable,2),[],'all')*1.1),-1);
o = semigeographicorder;
fignum = ["c.","a.","b."];
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
    c.TickLabels = aolabels;
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

