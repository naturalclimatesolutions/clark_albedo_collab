% Additional Data for 0.005 Maps

% Calculate Land Areas
% ********************
smtilesize = blocksize005 / 5;
ntilex = blocksize005/smtilesize;
ntiley = blocksize005/smtilesize;
nsmtiles = ntilex * ntiley;

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    tic
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'landmask','northlat','southlat','eastlon','westlon')
    
    pixarea = NaN(blocksize005,blocksize005);
    latser = northlat : - latlon005 : southlat;
    lonser = westlon : latlon005 : eastlon;
    
    % Because this function is slow, I'll only calculate areas where needed
    [a,b] = find(landmask);
    
    if numel(a) <= (smtilesize*smtilesize)
        for pt = 1 : numel(a)
            ii = a(pt);
            jj = b(pt);
            lat1 = latser(ii); lat2 = latser(ii+1);
            lon1 = lonser(jj); lon2 = lonser(jj+1);
            pixarea(ii,jj) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
        end
        save(subfilename,'pixarea','-append')
    else
        tiledlandmask = false(smtilesize,smtilesize,nsmtiles);
        tiledpixarea = nan(smtilesize,smtilesize,nsmtiles);
        tiledlats = zeros(smtilesize+1,nsmtiles);
        tiledlons = zeros(smtilesize+1,nsmtiles);
        tileland = true(nsmtiles,1);
        for tl = 1 : nsmtiles
            ii = smtilesi(tl,1) : smtilesi(tl,2);
            jj = smtilesj(tl,1) : smtilesj(tl,2);
            ttmask = landmask(ii,jj);
            if sum(ttmask,'all') == 0
                tileland(tl) = false;
            else
                tiledlandmask(:,:,tl) = ttmask;
                tiledpixarea(:,:,tl) = pixarea(ii,jj);
                i2 = smtilesi(tl,1) : smtilesi(tl,2)+1;
                j2 = smtilesj(tl,1) : smtilesj(tl,2)+1;
                tiledlats(:,tl) = latser(i2);
                tiledlons(:,tl) = lonser(j2);
            end
        end
        parfor tl = 1 : nsmtiles
            pres = tileland(tl);
            if ~pres, continue; end
            mask = tiledlandmask(:,:,tl);
            area = tiledpixarea(:,:,tl);
            lats = tiledlats(:,tl);
            lons = tiledlons(:,tl);
            [a,b] = find(mask);
            for pt = 1 : numel(a)
                ii = a(pt);
                jj = b(pt);
                lat1 = lats(ii); lat2 = lats(ii+1);
                lon1 = lons(jj); lon2 = lons(jj+1);
                area(ii,jj) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
            end
            tiledpixarea(:,:,tl) = area;
%             strcat("Done with tile #",num2str(tl))
        end
        for tl = 1 : nsmtiles
            ii = smtilesi(tl,1) : smtilesi(tl,2);
            jj = smtilesj(tl,1) : smtilesj(tl,2);
            pixarea(ii,jj) = tiledpixarea(:,:,tl);
        end
        save(subfilename,'pixarea','tileland','-append')
    end

    clear modisval northlat southlat westlon eastlon modislandmask pixarea latser lonser a b ...
        pt ii jj lat1 lat2 lon1 lon2 tiledpixarea tiledmodisval tiledlats tiledlons tiledlandmask ...
        tl ty tx i2 j2 

    strcat("Done with block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end


% Limit Bastin to areas that have 25% or higher opportunity and 25% or lower tree cover
% ************
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'bastin','bastinTP','treecover','treeselectedbastin')
    
    if exist('bastin','var') == 1
        bastinTP(bastinTP==bastinTPmiss) = nan;
        bastinCTC = bastinTP - bastin;
        if exist('treecover','var') == 1
            bastinCTC(bastinCTC<0) = treecover(bastinCTC<0);
            bastinCTC(isnan(bastinCTC)) = treecover(isnan(bastinCTC));
        else
            bastinCTC(bastinCTC<0) = 0;
        end

        selectedbastin = false(size(bastin));

        ii = bastin >= 25;
        jj = bastinCTC <= 25;
        kk = ii + jj == 2;
        selectedbastin(kk) = true;
        save(subfilename,'selectedbastin','bastinCTC','-append')
        min(bastin(kk))
    end
    if exist('treeselectedbastin','var') == 1
        treeselectedbastin = [];
        save(subfilename,'treeselectedbastin','-append')
    end
    clear subfilename bastin treecover bastinTP selectedbastin bastinCTC ii jj kk
    strcat("Done with selecting Bastin in block #",num2str(bb))
end
    

% Add Biomes and Continent arrays (to simplify statistics later)
% *******************************
bioperblock = false(nbblocks,nbiomes);
contperblock = false(nbblocks,numel(continents));
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,"ecoregion","region")

    biome = ones(blocksize005,blocksize005,'uint8') .* u8noval;
    thisblockeco = setdiff(unique(ecoregion),ecomis);
    for ee = 1 : numel(thisblockeco)
        eco = thisblockeco(ee);
        ecopix = ecoregion == eco;
        biome(ecopix) = biomeid(ecoid == eco);
    end
    thisblockbio = unique(biomeid(ismember(ecoid,thisblockeco)));
    bioperblock(bb,ismember(biomelist,thisblockbio)) = true;
    continent = ones(blocksize005,blocksize005,'uint8') .* u8noval;
    thisblockreglist = setdiff(unique(region),regsmis);
    for rr = 1 : numel(thisblockreglist)
        reg = thisblockreglist(rr);
        regpix = region == reg;
        continent(regpix) = contregid(reg);
    end
    thisblockcont = unique(contregid(thisblockreglist));
    contperblock(bb,thisblockcont) = true;

    save(subfilename,"biome","continent",'-append')
    strcat("Done with creating biome/contient array in block #",num2str(bb))

    clear biome continent ecoregion region thisblockreglist thisblockeco ee eco ecopix rr reg regpix
end

nodatabiome = u8noval;
nodatacontinent = u8noval;

save(sparameterfile,"nodatacontinent","nodatabiome","bioperblock","contperblock",'-append')


% Add combined mask (Griscom,Walker, Bastin)
% ******************************************
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,"selectedbastin","walkeroppcat","griscom")

    wlkm = ismember(walkeroppcat,Wsc);
    grim = griscom == 1;
    combinedopp = wlkm + grim + selectedbastin > 0;

    combinedoppclass = zeros(blocksize005,blocksize005,'uint8');
    combinedoppclass(wlkm) = combinedoppclass(wlkm) + 100;
    combinedoppclass(grim) = combinedoppclass(grim) + 10;
    combinedoppclass(selectedbastin) = combinedoppclass(selectedbastin) + 1;

    save(subfilename,"combinedopp","combinedoppclass",'-append')
    strcat("Done with creating combined opportunity masks in block #",num2str(bb))

    clear selectedbastin walkeroppcat griscom wlkm grim combinedopp combinedoppclass
end





