% Generate Most Likely Land cover

% Given that later I need them at resolution of 0.005, I will create it at that resolution.
% Also, now that I have two years of data, I will look through both to find the most likely, with
% priority to latest when there are two land cover of interest in the same spot.

% Note: my fist try was a disaster, as random (probably misclassified) pixels would translate over
%       large areas, so now I included a pixel count and use that, as well as a cutoff percent land
%       area to eliminate some land covers.


% Create indexes to simplify code
% *******************************
mlnames = ["MLOL","MLF","MLWSA"];
mllongnames = ["Most Likely Open Lands";"Most Likely Forests (exclusively)";...
    "Most Likely Forests and Woody Savannas";"Most Likely Forests"];
mlindnames = ["openlandind","forestind","forsavind"];
mllccomb = false(numel(mlnames),1); mllccomb(strcmp(mlnames,"MLWSA")) = true;
mleconames = strcat(mlnames,"perecoid");
mlbionames = strcat(mlnames,"perbiome");
ecoorder = sortrows(cat(2,(1:numel(ecoindex))',ecoid),2);
save(parameterfile,'mlnames','mllongnames','mlindnames','ecoorder','-append')

modisinorder = modisvars;
if strcmp(landcoverpriority,"newest")
    modisinorder = flip(modisinorder);
end    


% Extract most likely land cover per ecoregion (from land area statistics)
% ********************************************
% (Because it's easier to translate later, sort ecoregion by ID also)
EcoregionTotalLandArea = squeeze(sum(EcoLandcoverAreas(:,landind,:),2));
BiomesClimatesTotalLancoverAreas = sum(BiomesClimatesLancoverAreas(:,:,:,landind,:),4);
EcoregionTotalPixCount = squeeze(sum(EcoLandCoverPixelCounts(:,landind,2),2));
for ll = 1 : numel(mlnames)
    data = zeros(numel(ecoindex),1);
    eval(strcat("dataindex = ",mlindnames(ll),";"))
    for ee = 1 : numel(ecoindex)
        areas = squeeze(EcoLandcoverAreas(ee,dataindex,:));
        pixs = squeeze(EcoLandCoverPixelCounts(ee,dataindex,:));
        if mllccomb(ll) == true
            forest = areas(ismember(forsavind,forestind),:);
            sav = areas(ismember(forsavind,allsavind),:);
            csav = [sum(sav);zeros(numel(allsavind)-1,numel(modisvars))];
            areas = cat(1,forest,csav);
            fpix = pixs(ismember(forsavind,forestind),:);
            spix = pixs(ismember(forsavind,allsavind),:);
            cspix = [sum(spix);zeros(numel(allsavind)-1,numel(modisvars))];
            pixs = cat(1,fpix,cspix);
            clear forest sav csva fpix spix cspix
        end
        if ismember(biomeid(ee),tundradesertsindex)
            validareas = (areas ./ repmat(EcoregionTotalLandArea(ee,:),[numel(dataindex),1]))*100 ...
                >= prevalentlandcoverthreshold;
        else
            validareas = pixs >=minecopixels;
        end
        areas(logical(1-validareas)) = 0;
        maxarea = max(areas,[],'all');
        if maxarea == 0, continue; end
        [lctyp,yy] = find(areas == maxarea);
        val = unique(lctyp);
        if numel(val) > 1 && numel(unique(yy)) > 1
            if strcmp(mlnames(ll),"MLWSA")
                a = find(yy==min(yy));
                val = lctyp(a);
                if numel(val) > 1, error("check this"); end
            else
                a = find(yy==max(yy));
                val = lctyp(a);
                if numel(val) > 1, error("check this"); end
            end
        elseif numel(val) > 1 && numel(yy) == 1
            if strcmp(mlnames(ll),"MLWSA") && ismember(dataindex(lctyp),wsaind)
                dum = dataindex(lctyp);
                dum2 = setdiff(dum,wsaind);
                val = lctyp(dum==dum2);
                if numel(val) > 1, error("check this"); end
            else
                error("check this")
            end
        end
        data(ee) = dataindex(val);    
    end
    dum = data(ecoorder(:,1));
    dataperid = dum(2:numel(ecoid));   % first eco ID is zero, so I have to delete it!
    eval(strcat(mlnames(ll),"pereco = data;"))
    eval(strcat(mlnames(ll),"perecoid = dataperid;"))
    eval(strcat('mleconames(ll) = "',mlnames(ll),'perecoid";'))
    clear data areas maxareas ee
    
    data = zeros(numel(biomelist),numel(regionnames),numel(kglegnames));
    for bb = 1 : numel(biomelist)
        for rr = 1 : numel(regionnames)
            for cc = 1 : numel(kglegnames)
                areas = squeeze(BiomesClimatesLancoverAreas(bb,rr,cc,dataindex,:));
                pixs = squeeze(BiomesClimatesLandCoverPixelcounts(bb,rr,cc,dataindex,:));
                if mllccomb(ll) == true
                    forest = areas(ismember(forsavind,forestind),:);
                    sav = areas(ismember(forsavind,allsavind),:);
                    csav = [sum(sav);zeros(numel(allsavind)-1,numel(modisvars))];
                    areas = cat(1,forest,csav);
                    fpix = pixs(ismember(forsavind,forestind),:);
                    spix = pixs(ismember(forsavind,allsavind),:);
                    cspix = [sum(spix);zeros(numel(allsavind)-1,numel(modisvars))];
                    pixs = cat(1,fpix,cspix);
                    clear forest sav csva fpix spix cspix
                end
                if ismember(biomelist(bb),tundradesertsindex)
                    validareas = (areas ./ ...
                        squeeze(repmat(BiomesClimatesTotalLancoverAreas(bb,rr,cc,:,:),...
                        [1,1,1,numel(dataindex),1])))*100 >= prevalentlandcoverthreshold;
                else
                    validareas = pixs >=minbiomepixels;
                end
                areas(logical(1-validareas)) = 0;
                maxarea = max(areas,[],'all');
                if maxarea == 0, continue; end
                [lctyp,yy] = find(areas == maxarea);
                val = unique(lctyp);
                if numel(val) > 1 && numel(unique(yy)) > 1
                    if strcmp(mlnames(ll),"MLWSA")
                        a = find(yy==min(yy));
                        val = lctyp(a);
                        if numel(val) > 1, error("check this"); end
                    else
                        a = find(yy==max(yy));
                        val = lctyp(a);
                        if numel(val) > 1, error("check this"); end
                    end
                elseif numel(val) > 1 && numel(yy) == 1
                    if strcmp(mlnames(ll),"MLWSA") && ismember(dataindex(lctyp),wsaind)
                        dum = dataindex(lctyp);
                        dum2 = setdiff(dum,wsaind);
                        val = lctyp(dum==dum2);
                        if numel(val) > 1, error("check this"); end
                    else
                        error("check this")
                    end
                end               
                data(bb,rr,cc) = dataindex(val);
            end
        end
    end
    eval(strcat(mlnames(ll),"perbiome = data;"))
    eval(strcat('mlbionames(ll) = "',mlnames(ll),'perbiome";'))
    clear data areas maxareas ee dataindex
    
end
save(EcoStatfname,mleconames{:},mlbionames{:},'mleconames','mlbionames','ecoorder',...
    'mlnames','mlindnames','-append')


% Create a suite of different resolution gridded maps
% ***************************************************
gridlist = strings(numel(gridres),numel(mlnames),1);
for rr = 1 : numel(gridres)
    resolution = gridres(rr);
    rn = replace(num2str(resolution),".","_");
    for ll = 1 : numel(mlnames)
        eval(strcat('gridlist(rr,ll) = "',mlnames(ll),"GR",rn,'";'))
    end
end
        
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    tic
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    mlffname =  strcat(regoutputfiles,"MostLikely005_",num2str(bb),".mat");
    load(subfilename,modisvars{:},'ecoregion','pixarea','landmask')

    for yy = 1 : numel(modisvars)
        eval(strcat("landcover.",modisinorder(yy)," = ",modisinorder(yy),";"))
    end

    for ll = 1 : numel(mlnames)
        eval(strcat("landcoverindex = ",mlindnames(ll),";"))
        eval(strcat("mostlikelypereco = ",mleconames(ll),";"))
        combinesavannaflag = mllccomb(ll);
        prevalentlc = ones(blocksize005,blocksize005,numel(gridres),'uint8') .* u8noval;

        for rr = 1 : numel(gridres)     % coded like this to use in parfor
            resolution = gridres(rr);
            map = prevalentlandcovergrid(resolution,landcover,ecoregion,pixarea,landmask,...
                landcoverindex,mostlikelypereco,prevalentlandcoverthreshold,combinesavannaflag,...
                latlon005,ecomis);
            prevalentlc(:,:,rr) = map;
        end

        for rr = 1 : numel(gridres)
            res = gridres(rr);
            eval(strcat(mlnames(ll),"GR",replace(num2str(res),".","_")," = prevalentlc(:,:,rr);"))
        end
    end

    save(mlffname,gridlist{:})

    strcat("Done with Gridded Maps calculation in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end
save(sparameterfile,"gridlist",'-append')


% Create Most Likely Maps
% ***********************
potentialfix = zeros(500,6); ct = 1;
analysismaps = ["mask";"ecolevel";"biolevel";"onlylocalbiome";"allclimbiome";"otherrules","alevels"];
intmaps = strcat(repmat(mlnames,[numel(analysismaps),1]),repmat(analysismaps,[1,numel(mlnames),1]));
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    tic
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,modisvars{:},'ecoregion','landmask','region','koppengeiger')
    mlffname =  strcat(regoutputfiles,"MostLikely005_",num2str(bb),".mat");
    load(mlffname,gridlist{:})

    blockeco = ecoperblock(bb,:); blockeco = blockeco(isfinite(blockeco));

    for ll = 1 : numel(mlnames)
        eval(strcat("landcoverindex = ",mlindnames(ll),";"))
        data = ones(blocksize005,blocksize005,'uint8') .* u8noval;
        data(landmask) = 0;
        mask = false(blocksize005,blocksize005);
        leveldata = data;
        ii = landmask;
        % First fill all pixels where desired land cover exists
        for yy = 1 : numel(modisinorder)
            eval(strcat("modis =",modisinorder(yy),";"))
            if yy == 1, modisval = modis; end
            if mllccomb(ll) == true
                modis(modis==9) = 8; %#ok<SAGROW>
            end
            jj = ismember(modis,landcoverindex);
            kk = ii+jj == 2;
            data(kk) = modis(kk);
            ii = data == 0; clear kk
        end
        % Loop through all grids until desired land cover is found
        for rr = 1 : numel(gridres)
            eval(strcat("griddata = ",gridlist(rr,ll),";"))
            jj = ismember(griddata,landcoverindex);
            kk = ii + jj == 2;
            data(kk) = griddata(kk);
            leveldata(kk) = rr;
            ii = data == 0; clear jj kk
        end
        % For the pixels that still don't have a value, get the ecoregion most prevalent
        if sum(ii,'all') > 0
            eval(strcat("mlpeco = ",mleconames(ll),";"))
            for ee = 1 : numel(blockeco)
                thiseco = blockeco(ee);
                if blockeco(ee) == 0, continue; end
                jj = ecoregion == thiseco;
                kk = ii + jj == 2;
                data(kk) = mlpeco(thiseco);
                leveldata(kk) = numel(gridres)+10;
                ii = data == 0; clear jj kk
            end
        end
        % For the pixels that still don't have a value, get the biome/climate/region value
        ecodata = data;
        biodata = data;
        justbio = ones(blocksize005,blocksize005,'uint8') .* u8noval;
        justbio(landmask) = 0;
        justbiomore = justbio;
        otherrules = justbio;
        mask(ii) = true;   % save pixels that don't have values at this point
        [a,b] = find(ii);
        if numel(a) > 0
            eval(strcat("mlpbio = ",mlbionames(ll),";"))
            for pt = 1 : numel(a)
                ii = a(pt);
                jj = b(pt);
                eco = ecoregion(ii,jj);
                bio = biomeid(ecoid==eco);
                clim = koppengeiger(ii,jj);
                reg = region(ii,jj);
                if ismember(eco,[0,ecomis]) && ismember(clim,polardesertindex)
                    data(ii,jj) = 0;
                elseif ~isnan(bio)
                    if ismember(bio,tundradesertsindex) && ismember(clim,polardesertindex)
                        data(ii,jj) = 0;
                    elseif clim ~= misskoppengeiger && reg ~= regsmis
                        val = mlpbio(bio,reg,clim);
                        biodata(ii,jj) = val;
                        justbio(ii,jj) = val;
                        leveldata(ii,jj) = numel(gridres)+20;
                        if val == 0
                            dum = unique(mlpbio(bio,:,clim));
                            dum2 = dum(dum>0);
                            if numel(dum2) == 1
                                val = dum2;
                                justbiomore(ii,jj) = val;
                                leveldata(ii,jj) = numel(gridres)+30;
                            elseif numel(dum2) > 1
                                aa = zeros(numel(dum2),1);
                                for z = 1 : numel(dum2)
                                    sr = squeeze(mlpbio(bio,:,clim)) == dum2(z);
                                    aa(z) = sum(BiomesClimatesLancoverAreas...
                                        (bio,sr,clim,landcoverindex,2),'all');
                                end
                                m = aa == max(aa); 
                                val = dum2(aa == max(aa));
                                justbiomore(ii,jj) = val;
                                leveldata(ii,jj) = numel(gridres)+30;
                            elseif ismember(bio,tundradesertsindex) || ismember(clim,polardesertindex)
                                val = 0;
                            elseif contains(mlnames(ll),"WSA") && ...
                                    contains(biomenames(bio),"savanna",'IgnoreCase',true)
                                val = wsaind;
                                otherrules(ii,jj) = val;
                                leveldata(ii,jj) = numel(gridres)+40;
                            elseif contains(mlnames(ll),"WSA") && ...
                                    contains(econames(ecoid==eco),"savanna",'IgnoreCase',true)
                                val = wsaind;
                                otherrules(ii,jj) = val;
                                leveldata(ii,jj) = numel(gridres)+40;
                            elseif bio==2 && contains(kglegnames(clim),"Temperate","IgnoreCase",true)
                                val = LC005ind(contains(IGBPBiomes,"deciduous broadleaf",...
                                    "IgnoreCase",true));
                                otherrules(ii,jj) = val;
                                leveldata(ii,jj) = numel(gridres)+40;
                            else
                                pfb = potentialfix(:,1) == bb;
                                pfl = potentialfix(:,2) == ll;
                                pfo = potentialfix(:,3) == bio;
                                pfr = potentialfix(:,4) == double(reg);
                                pfc = potentialfix(:,5) == double(clim);
                                pfe = potentialfix(:,6) == double(eco);
                                if sum((pfb+pfl+pfo+pfr+pfc+pfe) == 6) == 0
                                    warning(strcat("Have a value of zero for ",mllongnames(ll),...
                                        ' in "',biomenames(bio),'" & "',kglegnames(clim),'"',...
                                        " found in block #",num2str(bb)))
                                    potentialfix(ct,:) = [bb,ll,bio,double(reg),double(clim),double(eco)];
                                    ct = ct+1;
                                end
                                val = 0;
                            end
                        end
                        data(ii,jj) = val;
                    end
                end
            end
        end
        eval(strcat(mlnames(ll)," = data;"))
        eval(strcat(mlnames(ll),"mask = mask;"))
        eval(strcat(mlnames(ll),"ecolevel = ecodata;"))
        eval(strcat(mlnames(ll),"biolevel = biodata;"))
        eval(strcat(mlnames(ll),"onlylocalbiome = justbio;"))
        eval(strcat(mlnames(ll),"allclimbiome = justbiomore;"))
        eval(strcat(mlnames(ll),"otherrules = otherrules;"))
        eval(strcat(mlnames(ll),"alevels = leveldata;"))

        clear landcoverindex data kk ii jj a b rr ee pt eco bio clim reg mlpeco mlpbio
    end
    

    % Add the combined most likely forest (and a maps of savanna locations)
    forestsavannamap = ones(blocksize005,blocksize005,'uint8') .* u8noval;
    forestsavannamap(landmask) = 1;
    MLFC = MLF;
    MLFCmask = MLFmask;
    MLFCecolevel = MLFecolevel;
    MLFCbiolevel = MLFbiolevel;
    MLFConlylocalbiome = MLFonlylocalbiome;
    MLFCallclimbiome = MLFallclimbiome;
    MLFCotherrules = MLFotherrules;
    MLFCalevels = MLFalevels;
    thisblocksavanna = setdiff(blockeco,forestedecoregions);
    for ee = 1 : numel(thisblocksavanna)
        kk = ecoregion == thisblocksavanna(ee);
        MLFC(kk) = MLWSA(kk);
        MLFCmask(kk) = MLWSAmask(kk);
        MLFCecolevel(kk) = MLWSAecolevel(kk);
        MLFCbiolevel(kk) = MLWSAbiolevel(kk);
        MLFConlylocalbiome(kk) = MLWSAonlylocalbiome(kk);
        MLFCallclimbiome(kk) = MLWSAallclimbiome(kk);
        MLFCotherrules(kk) = MLWSAotherrules(kk);
        MLFCalevels(kk) = MLWSAalevels(kk);
        forestsavannamap(kk) = 2;
    end

    % Add a map of forest/savanna

    addmaps = strcat("MLFC",analysismaps);
    
    save(mlffname,mlnames{:},intmaps{:},'MLFC',addmaps{:},'forestsavannamap','-append')

    clear(intmaps{:},mlnames{:},"MLFC",addmaps{:},"forestsavannamap")

    strcat("Done with Most Likely Map Calculation in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end

potentialfix = potentialfix(1:ct-1,:);

mlnames = [mlnames,"MLFC"];
mlindnames = [mlindnames,"forwsaind"];
save(sparameterfile,'mlnames','mlindnames','-append')


% For answer to reviewers (histogram of most-likely origin)
levlist = [1:numel(gridres),numel(gridres)+10,numel(gridres)+20,numel(gridres)+30,numel(gridres)+40];
MLOLdt = zeros(nbiomes,numel(gridres)+4);
MLFCdt = zeros(nbiomes,numel(gridres)+4);
for bio = 1 : nbiomes
    blocks = find(bioperblock(:,bio));
    histpblol = zeros(numel(blocks),numel(gridres)+4);
    histpblfs = zeros(numel(blocks),numel(gridres)+4);
    for bid = 1 : numel(blocks)
        bb = blocks(bid);
        load(strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat"),'biomes')
        load(strcat(regoutputfiles,"MostLikely005_",num2str(bb),".mat"),...
            'MLOLalevels','MLFCalevels')

        biomepix = biomes == biolist(bio);
        histpblol(bid,:) = histcounts(MLOLalevels(biomepix),[levlist,100]);
        histpblfs(bid,:) = histcounts(MLFCalevels(biomepix),[levlist,100]);
        clear biomes MLOLalevels MLFCalevels bb biomepix
    end
    MLOLdt(bio,:) = sum(histpblol);
    MLFCdt(bio,:) = sum(histpblfs);
    clear blocks histpblfs histpblol bid
end





    


    
            
    
        
    

    
