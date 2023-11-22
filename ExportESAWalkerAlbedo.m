% Export ESA-Truncated Walker Albedo-Combination
% **********************************************

trunc = ["walkertot80","walkertot85","walkertot90","walkertot95"];
varlist = [trunc,"WalkertotCO2"];
thres = ["wlk80";"wlk85";"wlk90";"wlk95"];

% Calculate Area statistics
% *************************
cvals = -2000:2000; ncvals = numel(cvals);
truncareacountarraylist = permute(strcat(varlist,"_AreaCount"),[2,1]);
blockarrays = strcat("Block_",truncareacountarraylist);
pfsize = [blocksize005,blocksize005,numel(truncareacountarraylist)];

for aa = 1 : numel(truncareacountarraylist)
    eval(strcat(truncareacountarraylist(aa)," = zeros(ncvals,1);"))
end

for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    statfilename = strcat(regoutputfiles,"ROstats005_",num2str(bb),".mat");
    load(subfilename,"pixarea","landmask")
    load(radforfname,varlist{:})

    data = nan(pfsize,'single');
    intarray = zeros(ncvals,numel(varlist));

    for aa = 1 : numel(varlist)
        eval(strcat("dum = ",varlist(aa),";"));
        data(:,:,aa) = round(dum);
    end
    
    parfor aa = 1 : numel(varlist)
        ardata = data(:,:,aa);
        ardata(landmask==0) = nan;
        ararray = intarray(:,aa);
        arpixarea = pixarea;
        datalist = unique(ardata(isfinite(ardata)));
        for cc = 1 : numel(datalist)
            cid = cvals == datalist(cc);
            pix = ardata == datalist(cc);
            ararray(cid) = sum(arpixarea(pix),'all','omitnan');
        end
        intarray(:,aa) = ararray;
    end
    
    for aa = 1 : numel(truncareacountarraylist)
        aname = blockarrays(aa);
        array = intarray(:,aa);
        eval(strcat(aname," = array;"))
        eval(strcat(truncareacountarraylist(aa)," = ",truncareacountarraylist(aa)," + array;"))
        clear aname array
    end
    clear aa masks data intarray
    
    % test that I get the same results as before
    load(statfilename,"Block_CarbonOnly_AreaCount")
    if exist("Block_CarbonOnly_AreaCount",'var')
        if ~isempty(find(abs(Block_CarbonOnly_AreaCount-Block_WalkertotCO2_AreaCount)>10-4,1))
            error(strcat("check for pixel counts in block #",num2str(bb)))
        end
    else
        Block_CarbonOnly_AreaCount = Block_WalkertotCO2_AreaCount;
        save(statfilename,"Block_CarbonOnly_AreaCount",'-append')
    end

    save(statfilename,blockarrays{1:numel(trunc)},'-append')
    clear(blockarrays{:},varlist{:},"pixarea","Block_CarbonOnly_AreaCount")

    strcat("Done with adding areas for truncated map stats in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
end

truncstatsarraylist = permute(strcat(truncareacountarraylist,"Stats"),[2,1]);
for aa = 1 : numel(truncareacountarraylist)
    eval(strcat("data = ",truncareacountarraylist(aa),";"))
    vals = cvals;
    if contains(truncareacountarraylist(aa),"AlbedoOffset"), vals = aovals; end
    statsarray = zeros(numel(mappct),1);
    totarea = sum(data);
    cumdata = cumsum(data);
    areastval = totarea .* mappct / 100;
    for ss = 1 : numel(mappct)
        ii = find(cumdata>=areastval(ss),1);
        statsarray(ss) = vals(ii); clear ii
    end
    eval(strcat(truncstatsarraylist(aa)," = statsarray;"))
    clear data vals statsarray totarea cumdata areastval ss
end

save(ReforOppfname, truncstatsarraylist{:},truncareacountarraylist{:},'-append')
save(sparameterfile,"truncareacountarraylist","truncstatsarraylist",'-append')


% Print truncated maps
% ********************
indexfileroot = strcat(regoutputfiles,"ROinputs_");
fileroot = strcat(regoutputfiles,"RadForcing005_");
inputdataformat = "single";
lowerresolution = 0.1;
method = "mean";
missinglandblock = 51;

latlr = latlim005(2) - lowerresolution/2 : -lowerresolution : latlim005(1) + lowerresolution/2;
lonlr = lonlim005(1) + lowerresolution/2 : lowerresolution : lonlim005(2) - lowerresolution/2;
[lons,lats] = meshgrid(lonlr,latlr);

thisfiguretype = "three panels";

for ll = 1 : numel(trunc)
    input = trunc(ll);
    [lowresmap,~] = printmapprep(input,fileroot,preswlk,latlon005,inputdataformat,...
        rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot);
    
    lrvarname = strcat(input,"lowres");
    eval(strcat(lrvarname," = lowresmap;"));
    
    save(ReforOppfname,lrvarname{:},'-append')
    
    axislegend = 'Mg CO_2e ha^-^1';
    
    eval(strcat("prctvals =",truncstatsarraylist(ll),";"))
    
    % Pring jpeg figure
    figurename = strcat(figuredir,input);
    figure(ll); clf
    h = albedofigure(lowresmap,seqcat,lats,lons,figurename,thisfiguretype,...
    co2colorbar,prctvals,coastlat,coastlon,axislegend,false);
    
    clear lowresmap highresmap
end



% Alternative maps (only look into 85% for now)
% ****************
intype = ["NCI","AO"];
vars.base = "NCIbase";
vars.modified = "NCIwlk85";
copmmet = ["difference","percent"];
thisfiguretype = "two panels";
for tt = 1 : numel(intype)
    vars.base = strcat(intype(tt),"base");
    vars.modified = strcat(intype(tt,"wlk85"));
    for mm = 1 : numel(copmmet)
        vars.method = copmmet(mm);
        [lowresmap,~] = printmapprep(vars,fileroot,preswlk,latlon005,inputdataformat,...
            rastersize005,blocksize005,lowerresolution,method,missinglandblock,indexfileroot);

        % I should not have to do this, but I don't have time now to seek out the few pixels that
        % have minimal positiva values
        lowresmap(lowresmap > 0) = 0;

        lrvarname = strcat("NCIWlk85",copmmet(mm),"lowres");
        eval(strcat(lrvarname," = lowresmap;"));
        save(ReforOppfname,lrvarname,'-append')

        if strcmp(copmmet,"difference")
            if strcmp(intype,"NCI")
                thesecat = seqcat;
                thiscolorbar = co2colorbar;
                axislegend = 'Mg CO_2e ha^-^1';
            else
                thesecat = offsetcat;
                thiscolorbar = aocolorbar;
                axislengend = "%";
            end
        else
            axislegend = 'Percent Differences [%]';
            thesecat = [-1000,-50,-25,-10,-5,0];
            thiscolorbar = flip([254,235,226;251,180,185;247,104,161;197,27,138;122,1,119]/255);
%             thiscolorbar = flip([254,235,226;252,197,192;250,159,181;...
%                 247,104,161;221,52,151;174,1,126;122,1,119]/255);
        end

        figurename = strcat(figuredir,"Walker85",intype(tt),"Diff",copmmet(mm));
        figure(20+tt+mm); clf
        h = albedofigure(lowresmap,thesecat,lats,lons,figurename,thisfiguretype,thiscolorbar,...
            [],coastlat,coastlon,axislegend,false);
    end
end


% Export Truncated stats in opportunity tables
% ********************************************
mapcats = ["AO","NCI"];
nopp = numel(allreforestationopp);
inputs = reshape(strcat(repmat(mapcats,[numel(thres),1]),...
    repmat(thres,[1,numel(mapcats)])),[numel(thres)*numel(mapcats),1]);
dum = strcat("AO",thres); % Now only calculate stats by AO bins ...
truncareatables = strcat("AreasbyBiome_",dum);
truncco2tables = strcat("TotalCO2byBiome_",dum);
truncseqtables = strcat("JustCarbonbyBiome_",dum);

aoarray = zeros(nbiomes,nocat,nopp);
for tt = 1 : numel(truncareatables)    
    eval(strcat(truncareatables(tt)," = aoarray;"))
    eval(strcat(truncco2tables(tt)," = aoarray;"))
    eval(strcat(truncseqtables(tt)," = aoarray;"))
end
Btareatables = strcat("Block",truncareatables);
Btco2tables = replace(Btareatables,"Areas","TotalCO2");
Btseqtables = replace(Btareatables,"Areas","JustCarbon");
negco2 = offsetcat(1:nocat)>=100;

% Calculate Statistics
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    statfilename = strcat(regoutputfiles,"ROstats005_",num2str(bb),".mat");
    load(radforfname,'RFmed005')

    if sum(isnan(RFmed005),'all') == blocksize005^2, clear('RFmed005'), continue, end
    
    load(subfilename,"biome","continent","pixarea","selectedbastin","walkeroppcat","griscom","landmask")
    load(radforfname,inputs{:},trunc{:},"AObase")
    
    areaarray = repmat(aoarray,[1,1,1,numel(thres)]);
    co2array = areaarray;
    seqarray = areaarray;
    allnci = nan(blocksize005,blocksize005,numel(thres));
    allao = nan(blocksize005,blocksize005,numel(thres));
    allcarb = nan(blocksize005,blocksize005,numel(thres));

    for ss = 1 : numel(thres)
        eval(strcat("allcarb(:,:,ss) = ",trunc(ss),";"));
        eval(strcat("allnci(:,:,ss) = NCI",thres(ss),";"))
        eval(strcat("allao(:,:,ss) = AO",thres(ss),";"))
    end

    parfor ss = 1 : numel(thres)
        nci = allnci(:,:,ss);
        ao = allao(:,:,ss);
        carb = allcarb(:,:,ss);
        sarea = areaarray(:,:,:,ss);
        mitigation = co2array(:,:,:,ss);
        noalbcar = seqarray(:,:,:,ss);

        for oo = 1 : nopp
            refopp = allreforestationopp(oo); %#ok<PFBNS> 
            switch refopp
                case "Globally"
                    thismask = landmask;
                case "Walker"
                    thismask = ismember(walkeroppcat,Wsc);
                case "Griscom"
                    thismask = griscom == 1;
                case "Bastin"
                    thismask = selectedbastin;
                case "CombinedOpp"
                    opmask = combinedopp;
            end

            [area,co2,seq,missingareas,mismask] = binnedvalues(ao,carb,nci,...
                thismask,offsetcat,negco2,continent,biome,pixarea,continentlist,biomelist,...
                nodatacontinent,nodatabiome,landmask);
            sarea(:,:,oo) = area(:,:,ncont+1);
            mitigation(:,:,oo) = co2(:,:,ncont+1);
            noalbcar(:,:,oo) = seq(:,:,ncont+1);
        end

        areaarray(:,:,:,ss) = sarea;
        co2array(:,:,:,ss) = mitigation;
        seqarray(:,:,:,ss) = noalbcar;


    end

    for tt = 1 : numel(thres)
        eval(strcat("Block",truncareatables(tt)," = areaarray(:,:,:,tt);"))
        eval(strcat("Block",truncco2tables(tt)," = co2array(:,:,:,tt);"))
        eval(strcat("Block",truncseqtables(tt)," = seqarray(:,:,:,tt);"))
        eval(strcat(truncareatables(tt)," = ",truncareatables(tt)," + Block",truncareatables(tt),";"))
        eval(strcat(truncco2tables(tt)," = ",truncco2tables(tt)," + Block",truncco2tables(tt),";"))
        eval(strcat(truncseqtables(tt)," = ",truncseqtables(tt)," + Block",truncseqtables(tt),";"))
    end

    load(statfilename,"BlockAreasbyBiome_AObase")
    bioareas = squeeze(sum(BlockAreasbyBiome_AObase(:,:,:,ncont+1),2));
    [nba,noa] = size(bioareas);
    for tt = 1 : numel(thres)
        trareas = squeeze(sum(areaarray(:,:,:,tt),2));
        [nbt,nott] = size(trareas);
        if noa < nott
            ii = ~ismember(allreforestationopp,"CombinedOpp");
            trareas = trareas(:,ii);
        end
        if nba < nbt
            ii = ~ismember(biomenames,["N/A","Rock and Ice"]);
            trareas = trareas(ii,:);
        end
        if ~isempty(find(abs(bioareas-trareas)>10^-4,1))
            error(strcat("My sums are not the same as before in block #",num2str(bb)))
        end
    end
    clear BlockAreasbyBiome_AObase bioareas tt trareas

    save(statfilename,Btareatables{:},Btco2tables{:},Btseqtables{:},'-append')
    
    clear(inputs{:},trunc{:},"biome","pixarea","selectedbastin","walkeroppcat","griscom","landmask",'RFmed005')
    clear(Btareatables{:},Btco2tables{:},Btseqtables{:})
    
    strcat("Done with sensitivity AO/NCI statistics block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end

save(ReforOppfname, truncareatables{:},truncco2tables{:},truncseqtables{:},"MissRFareas",'-append')
save(sparameterfile,"sensareatables","sensco2tables","sensseqtables","bdnames",...
    "allreforestationopp",'-append')


% Export tables
% *************
below50 = offsetcat(1:nocat)<aothreshold;
above75 = offsetcat(2:nocat+1)>=75;
above100 = offsetcat(2:nocat+1)>=100;
semigeographicorder = [11,6,8,12,5,4,7,3,2,1,14,9,10,13];
opporder = ["Griscom","Bastin","Walker","TotalBiome"];
no = numel(opporder);
[~,noi] = ismember(opporder,allreforestationopp);
noi(noi==0) = 4;
ncont = numel(continents);
clm = ["Biomes","Area","Or. carbon","Trunc carbon","Carbon change",...
    "Or. Walker NCI","trunc NCI","NCI reduction",...
    "<Original 50%AO %Areas","trunc <50%AO %Areas","<50% Areas decrease",...
    "Original >100%AO %Areas","trunc >100%AO %Areas","Climate negative increase"];

clu = ["","[Mha]","[Pg CO2e]","[Pg CO2e]","[%]","[Pg CO2e]","[Pg CO2e]","[%]",...
    "[%]","[%]","[%]","[%]","[%]","[%]"];

worldvalues = zeros((nbiomes+1)*(no+1),numel(clm)-1);

abase = AreasbyBiome_AObase(semigeographicorder,:,noi,ncont+1);
cbase = TotalCO2byBiome_AObase(semigeographicorder,:,noi,ncont+1);
sbase = JustCarbonbyBiome_AObase(semigeographicorder,:,noi,ncont+1);

[~,~,tnop,~] = size(AreasbyBiome_AOwlk85);
if tnop == no+1
    ii = ~ismember(allreforestationopp,"CombinedOpp");
    rarea = AreasbyBiome_AOwlk85(:,:,ii);
    rco2 = TotalCO2byBiome_AOwlk85(:,:,ii);
    rseq = JustCarbonbyBiome_AOwlk85(:,:,ii);
end
trarea = rarea(semigeographicorder,:,noi);
trco2 = rco2(semigeographicorder,:,noi);
trseq = rseq(semigeographicorder,:,noi);

for bb = 1 : nbiomes
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,1) = sum(abase(bb,:,:),2);
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,2) = sum(sbase(bb,:,:),2);
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,3) = sum(trseq(bb,:,:),2);
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,4) = abs(sum(sbase(bb,:,:),2) - sum(trseq(bb,:,:),2)) ./sum(sbase(bb,:,:),2) .* 100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,5) = sum(cbase(bb,:,:),2);
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,6) = sum(trco2(bb,:,:),2);
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,7) = abs(sum(cbase(bb,:,:),2) - sum(trco2(bb,:,:),2)) ./sum(cbase(bb,:,:),2) .* 100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,8) = sum(abase(bb,below50,:),2) ./sum(abase(bb,:,:),2) .*100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,9) = sum(trarea(bb,below50,:),2) ./sum(trarea(bb,:,:),2) .*100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,10) = (sum(trarea(bb,below50,:),2) ./sum(trarea(bb,:,:),2) - sum(abase(bb,below50,:),2) ./sum(abase(bb,:,:),2)) .*100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,11) = sum(abase(bb,above100,:),2) ./sum(abase(bb,:,:),2) .*100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,12) = sum(trarea(bb,above100,:),2) ./sum(trarea(bb,:,:),2) .*100;
    worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,13) = (sum(trarea(bb,above100,:),2) ./sum(trarea(bb,:,:),2) - sum(abase(bb,above100,:),2) ./sum(abase(bb,:,:),2)) .*100;
end
fc = reshape(cat(1,biomenames(semigeographicorder)',repmat(opporder',[1,nbiomes])),[nbiomes*(no+1),1]);

worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,1) = sum(abase,1:2);
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,2) = sum(sbase,1:2);
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,3) = sum(trseq,1:2);
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,4) = abs(sum(sbase,1:2) - sum(trseq,1:2)) ./sum(sbase,1:2) .* 100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,5) = sum(cbase,1:2);
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,6) = sum(trco2,1:2);
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,7) = abs(sum(cbase,1:2) - sum(trco2,1:2)) ./sum(cbase,1:2) .* 100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,8) = sum(abase(:,below50,:),1:2) ./sum(abase,1:2) .*100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,9) = sum(trarea(:,below50,:),1:2) ./sum(trarea,1:2) .*100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,10) = (sum(trarea(:,below50,:),1:2) ./sum(trarea,1:2) - sum(abase(:,below50,:),1:2) ./sum(abase,1:2)) .*100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,11) = sum(abase(:,above100,:),1:2) ./sum(abase,1:2) .*100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,12) = sum(trarea(:,above100,:),1:2) ./sum(trarea,1:2) .*100;
worldvalues(nbiomes*(no+1)+2:nbiomes*(no+1)+no+1,13) = (sum(trarea(:,above100,:),1:2) ./sum(trarea,1:2) - sum(abase(:,above100,:),1:2) ./sum(abase,1:2)) .*100;

T = cell((nbiomes+1)*(no+1)+2,numel(clm));
T(1,:) = cellstr(clm);
T(2,:) = cellstr(clu);
T(3:nbiomes*(no+1)+2,1) = cellstr(fc);
T(nbiomes*(no+1)+3:(nbiomes+1)*(no+1)+2,1) = cellstr(cat(2,"Globally",opporder));
T(3:(nbiomes+1)*(no+1)+2,2:numel(clm)) = num2cell(worldvalues);

TT = cell2table(T);
writetable(TT,xcelfile,'Sheet','Truncated')






for ss = 1 : numel(thres)
	sensname = thres(ss);
    atbl = eval(strcat("AreasbyBiome_AO",sensname,"(:,:,noi);"));
    ctbl = eval(strcat("TotalCO2byBiome_AO",sensname,"(:,:,noi);"));
    stbl = eval(strcat("JustCarbonbyBiome_AO",sensname,"(:,:,noi);"));
    for bb = 1 : nbiomes
        bio = semigeographicorder(bb);
        tar = squeeze(sum(atbl(bio,:,:),2));
        mit = squeeze(sum(ctbl(bio,:,:),2));
        basetar = squeeze(sum(abase(bio,:,:),2));
        basemit = squeeze(sum(cbase(bio,:,:),2));
        posar = squeeze(sum(atbl(bio,below50,:),2));
        posmit = squeeze(sum(ctbl(bio,below50,:),2));
        posbasear = squeeze(sum(abase(bio,below50,:),2));
        posbasemit = squeeze(sum(cbase(bio,below50,:),2));
        
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,1,ss) = squeeze(sum(atbl(bio,:,:),2));
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,2,ss) = squeeze(sum(stbl(bio,:,:),2));
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,3,ss) = mit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,4,ss) = mit*co2scalar ./(tar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,5,ss) = basemit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,6,ss) = basemit*co2scalar ./(basetar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,7,ss) = basemit - mit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,8,ss) = (basemit - mit) ./ basemit .* 100;
        
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,9,ss) = posar./tar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,10,ss) = posmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,11,ss) = posmit*co2scalar ./(posar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,12,ss) = posbasear./basetar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,13,ss) = posbasemit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,14,ss) = posbasemit*co2scalar ./(posbasear*areascalar) * perareascalar;

        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,15,ss) = posbasemit - posmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,16,ss) = (posbasemit - posmit) ./ posbasemit .* 100;
        
        if ss == 1
            firstcolumn((bb-1)*(no+1)+1) = biomenames(bio);
            firstcolumn((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1) = opporder;
        end
    end
end

for ss = 1 : numel(trunc)
    sensname = trunc(ss);
    T = table(firstcolumn,worldvalues(:,1,ss),worldvalues(:,2,ss),worldvalues(:,3,ss),...
        worldvalues(:,4,ss),worldvalues(:,5,ss),worldvalues(:,6,ss),worldvalues(:,7,ss),...
        worldvalues(:,8,ss),worldvalues(:,9,ss),worldvalues(:,10,ss),worldvalues(:,11,ss),...
        worldvalues(:,12,ss),worldvalues(:,13,ss),worldvalues(:,14,ss),worldvalues(:,15,ss),...
        worldvalues(:,16,ss));
    T.Properties.VariableNames = clm;
    xcelfile = strcat(resfolder,"ReforestationOpportunityStatsbybiomes.xlsx");
    writetable(T,xcelfile,'Sheet',sensname)
end


% Combined Table
% **************
T = cell(nbiomes*(no+1)+2,numel(clm));
T(1,:) = cellstr(clm);
T(2,:) = cellstr(clu);
abase = AreasbyBiome_AObase(:,:,noi,ncont+1);
tar = squeeze(sum(abase,2)); posbar = squeeze(sum(abase(:,below50,:),2));
ppbar = posbar ./ tar .*100;
posiar = squeeze(sum(AreasbyBiome_AOminalb(:,below50,noi,ncont+1),2)); ppiar = posiar ./ tar .*100;
posaar = squeeze(sum(AreasbyBiome_AOmaxalb(:,below50,noi,ncont+1),2)); ppaar = posaar ./ tar .*100;
cbase = TotalCO2byBiome_AObase(:,:,noi,ncont+1);
bmit = squeeze(sum(cbase,2)); posbmit = squeeze(sum(cbase(:,below50,:),2));
bmitav = bmit.*co2scalar ./(tar.*areascalar) * perareascalar;
bposmitav = posbmit.*co2scalar ./(posbar.*areascalar) * perareascalar;
cmin = TotalCO2byBiome_AOminalb(:,:,noi,ncont+1);
imit = squeeze(sum(cmin,2)); posimit = squeeze(sum(cmin(:,below50,:),2));
imitav = imit.*co2scalar ./(tar.*areascalar) * perareascalar;
posimitav = posimit.*co2scalar ./(posiar.*areascalar) * perareascalar;
cmax = TotalCO2byBiome_AOmaxalb(:,:,noi,ncont+1);
amit = squeeze(sum(cmax,2)); posamit = squeeze(sum(cmax(:,below50,:),2));
amitav = amit.*co2scalar ./(tar.*areascalar) * perareascalar;
posamitav = posamit.*co2scalar ./(posaar.*areascalar) * perareascalar;
carb = squeeze(sum(JustCarbonbyBiome_AObase(:,:,noi,ncont+1),2));

for bb = 1 : nbiomes
    bio = semigeographicorder(bb);
    T((bb-1)*(no+1)+3,1) = cellstr(biomenames(bio));
    for oo = 1 : numel(opporder)
        T((bb-1)*(no+1)+3+oo,1) = cellstr(opporder(oo));
        T((bb-1)*(no+1)+3+oo,2) = num2cell(round(tar(bio,oo),2));
        T((bb-1)*(no+1)+3+oo,3) = num2cell(round(carb(bio,oo),2));
        T((bb-1)*(no+1)+3+oo,4) = cellstr(strcat(num2str(round(bmit(bio,oo),2),'%.2f')," [",...
            num2str(round(amit(bio,oo),2),'%.2f')," - ",num2str(round(imit(bio,oo),2),'%.2f'),"]"));
        T((bb-1)*(no+1)+3+oo,5) = cellstr(strcat(num2str(round(bmitav(bio,oo),2),'%.2f')," [",...
            num2str(round(amitav(bio,oo),2),'%.2f')," - ",num2str(round(imitav(bio,oo),2),'%.2f'),"]"));
        T((bb-1)*(no+1)+3+oo,6) = cellstr(strcat(num2str(round(ppbar(bio,oo)))," [",...
            num2str(round(ppaar(bio,oo)))," - ",num2str(round(ppiar(bio,oo))),"]"));
        T((bb-1)*(no+1)+3+oo,7) = cellstr(strcat(num2str(round(posbmit(bio,oo),2),'%.2f')," [",...
            num2str(round(posamit(bio,oo),2),'%.2f')," - ",num2str(round(posimit(bio,oo),2),'%.2f'),"]"));
        T((bb-1)*(no+1)+3+oo,8) = cellstr(strcat(num2str(round(bposmitav(bio,oo),2),'%.2f')," [",...
            num2str(round(posamitav(bio,oo),2),'%.2f')," - ",num2str(round(posimitav(bio,oo),2),'%.2f'),"]"));
    end
end

TT = cell2table(T);
writetable(TT,xcelfile,'Sheet','combined')

% Tables by opportunity
clm = ["Biomes","Area","Base NCI","Truncated NCI","% Reduction in NCI",...
    "%Area in AO lower than 50%","Reduction in Area of AO lower than 50%",...
    "%Area in AO higher than 75%","Increased areas in >75% AO",...
    "%Area in climate negative","Increased areas climate negative"];
clu = cat(2,["","[Mha]","[Pg CO2e]","[Pg CO2e]"],repmat("[%]",[1,7]));
atbl = AreasbyBiome_AOwlk85(semigeographicorder,:,noi);
ctbl = TotalCO2byBiome_AOwlk85(semigeographicorder,:,noi);
stbl = JustCarbonbyBiome_AOwlk85(semigeographicorder,:,noi);
abase = AreasbyBiome_AObase(semigeographicorder,:,noi,ncont+1);
cbase = TotalCO2byBiome_AObase(semigeographicorder,:,noi,ncont+1);
sbase = JustCarbonbyBiome_AObase(semigeographicorder,:,noi,ncont+1);

for ll = 1 : numel(opporder)
    sheetname = opporder(ll);
    worldvalues = zeros(nbiomes,numel(clm)-1);
    area = sum(atbl(:,:,ll),2);
    mit = sum(ctbl(:,:,ll),2);
    basemit = sum(cbase(:,:,ll),2);
    posarea = sum(atbl(:,below50,ll),2) ./ area .*100;
    posba = sum(abase(:,below50,ll),2) ./ area .*100;
    a75 = sum(atbl(:,above75,ll),2) ./ area .*100;
    ab75 = sum(abase(:,above75,ll),2) ./ area .*100;
    negarea = sum(atbl(:,above100,ll),2) ./ area .*100;
    negba = sum(abase(:,above100,ll),2) ./ area .*100;

    worldvalues(:,1) = area;
    worldvalues(:,2) = basemit;
    worldvalues(:,3) = mit;
    worldvalues(:,4) = (basemit - mit) ./basemit .*100;
    worldvalues(:,5) = posarea;
    worldvalues(:,6) = posba - posarea;
    worldvalues(:,7) = a75;
    worldvalues(:,8) = a75 - ab75;
    worldvalues(:,9) = negarea;
    worldvalues(:,10) = negarea - negba;
    
    T = cell(nbiomes+3,numel(clm));
    T(1,:) = cellstr(clm);
    T(2,:) = cellstr(clu);
    T(3:nbiomes+2,1) = cellstr(biomenames(semigeographicorder));
    T(3:nbiomes+2,2:numel(clm)) = num2cell(worldvalues);
    
    TT = cell2table(T);
    writetable(TT,xcelfile,'Sheet',sheetname)
        
end
    
    
    
    
    for bb = 1 : nbiomes
        bio = semigeographicorder(bb);
        worldvalues(bb+3,2) = sum(At(bio,:),2);
        worldvalues(bb+3,3) = sum(cbase(bio,:,ll),2);
        worldvalues(bb+3,4) = sum(Ct(bio,:),2);
        tar = squeeze(sum(At(bio,:,:),2));
        mit = squeeze(sum(Ct(bio,:,:),2));
        basetar = squeeze(sum(abase(bio,:,:),2));
        basemit = squeeze(sum(cbase(bio,:,:),2));
        posar = squeeze(sum(atbl(bio,below50,:),2));
        posmit = squeeze(sum(ctbl(bio,below50,:),2));
        posbasear = squeeze(sum(abase(bio,below50,:),2));
        posbasemit = squeeze(sum(cbase(bio,below50,:),2));
        
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,1,ss) = squeeze(sum(atbl(bio,:,:),2));
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,2,ss) = squeeze(sum(stbl(bio,:,:),2));
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,3,ss) = mit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,4,ss) = mit*co2scalar ./(tar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,5,ss) = basemit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,6,ss) = basemit*co2scalar ./(basetar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,7,ss) = basemit - mit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,8,ss) = (basemit - mit) ./ basemit .* 100;
        
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,9,ss) = posar./tar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,10,ss) = posmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,11,ss) = posmit*co2scalar ./(posar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,12,ss) = posbasear./basetar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,13,ss) = posbasemit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,14,ss) = posbasemit*co2scalar ./(posbasear*areascalar) * perareascalar;

        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,15,ss) = posbasemit - posmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,16,ss) = (posbasemit - posmit) ./ posbasemit .* 100;
        
        if ss == 1
            firstcolumn((bb-1)*(no+1)+1) = biomenames(bio);
            firstcolumn((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1) = opporder;
        end
    end

    
    


