% Calculate Opportunity Stats
% ***************************

% Prepare Arrays
% **************
nopp = numel(allreforestationopp);
Wsc = WlkClassID(contains(WlkClassDef,"R/"));
ncont = numel(continents);
continentlist = 1 : ncont;
sensareatables = strcat("AreasbyBiome_",AOarrays);
sensco2tables = strcat("TotalCO2byBiome_",AOarrays);
sensseqtables = strcat("JustCarbonbyBiome_",AOarrays);
MissRFareas = zeros(nbiomes+1,nopp);

aoarray = zeros(nbiomes,nocat,nopp,ncont+1);
for tt = 1 : numel(bdnames)    
    eval(strcat(sensareatables(tt)," = aoarray;"))
    eval(strcat(sensco2tables(tt)," = aoarray;"))
    eval(strcat(sensseqtables(tt)," = aoarray;"))
end

Bsareatables = strcat("Block",sensareatables);
Bsco2tables = replace(Bsareatables,"Areas","TotalCO2");
Bsseqtables = replace(Bsareatables,"Areas","JustCarbon");

negco2 = offsetcat(1:nocat)>=100;


% Calculate Statistics
% ********************
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    load(subfilename,"biome","continent","pixarea","selectedbastin",...
        "walkeroppcat","griscom","landmask","combinedopp")
    load(radforfname,"GWPmed005")
    statfilename = strcat(regoutputfiles,"ROstats005_",num2str(bb),".mat");

    if sum(isnan(GWPmed005),'all') == blocksize005^2
        mispix = zeros(nbiomes+1,nopp);
        for oo = 1 : nopp
            refopp = allreforestationopp(oo);
            switch refopp
                case "Globally"
                    opmask = landmask;
                case "Walker"
                    opmask = ismember(walkeroppcat,Wsc);
                case "Griscom"
                    opmask = griscom == 1;
                case "Bastin"
                    opmask = selectedbastin;
                case "CombinedOpp"
                    opmask = combinedopp;
            end
            if sum(opmask,'all') > 0
                thisbio = unique(biome(opmask));
                for k = 1 : numel(thisbio)
                    bio = thisbio(k); if bio==nodatabiome, bio = nbiomes+1; end
                    biopix = biome == bio;
                    kk = biopix + opmask == 2;
                    mispix(bio,oo) = sum(pixarea(kk),'all','omitnan') ./ areascalar;
                    if sum(isnan(pixarea(kk))) > 0
                        warning(strcat("See why I have missing pixel areas in block #",num2str(bb)))
                    end
                end
            end
        end
        if sum(mispix,'all') > 0
            MissRFareas = MissRFareas + mispix;
        end
        strcat("Done with sensitivity AO/NCI statistics block #",num2str(bb)," - no RF data")
        clear GWPmed005 combinedopp
        clear("biome","continent","pixarea","selectedbastin","walkeroppcat","griscom","landmask")
        continue
    end
    
    load(radforfname,AOarrays{:},NCIarrays{:},"WalkertotCO2")
    
    masks = false(blocksize005,blocksize005,nopp);
    for oo = 1 : nopp
        refopp = allreforestationopp(oo);
        switch refopp
            case "Globally"
                masks(:,:,oo) = landmask;
            case "Walker"
                masks(:,:,oo) = ismember(walkeroppcat,Wsc);
            case "Griscom"
                masks(:,:,oo) = griscom == 1;
            case "Bastin"
                masks(:,:,oo) = selectedbastin;
            case "CombinedOpp"
                masks(:,:,oo) = combinedopp;
        end
    end

    for ss = 1 : numel(bdnames)
        sensname = bdnames(ss);
        areaarray = permute(aoarray,[1,2,4,3]);
        co2array = areaarray;
        seqarray = areaarray;

        eval(strcat("nci = NCI",sensname,";"))
        eval(strcat("input = AO",sensname,";"))

        tmpmissarea = zeros(nbiomes+1,nopp);

        parfor oo = 1 : nopp
            thismask = masks(:,:,oo);
            
            [area,co2,seq,missingareas,mismask] = binnedvalues(input,WalkertotCO2,nci,...
                thismask,offsetcat,negco2,continent,biome,pixarea,continentlist,biomelist,...
                nodatacontinent,nodatabiome,landmask);
            areaarray(:,:,:,oo) = area;
            co2array(:,:,:,oo) = co2;
            seqarray(:,:,:,oo) = seq;
            tmpmissarea(:,oo) = sum(missingareas,2);
            
        end
        if ss == 1
            MissRFareas = MissRFareas + tmpmissarea;
        end
        eval(strcat("BlockAreasbyBiome_AO",sensname," = permute(areaarray,[1,2,4,3]);"))
        eval(strcat("BlockTotalCO2byBiome_AO",sensname," = permute(co2array,[1,2,4,3]);"))
        eval(strcat("BlockJustCarbonbyBiome_AO",sensname," = permute(seqarray,[1,2,4,3]);"))
    end

    for tt = 1 : numel(sensareatables)
        eval(strcat(sensareatables(tt)," = ",sensareatables(tt)," + Block",sensareatables(tt),";"))
        eval(strcat(sensco2tables(tt)," = ",sensco2tables(tt)," + Block",sensco2tables(tt),";"))
        eval(strcat(sensseqtables(tt)," = ",sensseqtables(tt)," + Block",sensseqtables(tt),";"))
    end

    save(statfilename,Bsareatables{:},Bsco2tables{:},Bsseqtables{:},'-append')
    
    clear(AOarrays{:},NCIarrays{:},'GWPmed005',"biome","continent",...
        "pixarea","WalkertotCO2","selectedbastin","walkeroppcat","griscom","landmask")
    clear(Bsareatables{:},Bsco2tables{:},Bsseqtables{:})
    
    strcat("Done with opportunity map AO/NCI statistics block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end


% Save arrays
% ***********
save(ReforOppfname, sensareatables{:},sensco2tables{:},sensseqtables{:},"MissRFareas",'-append')
save(sparameterfile,"sensareatables","sensco2tables","sensseqtables","bdnames",'-append')


% Create a sensitivity map (and calculate areas)
% ************************
sensarreas50 = zeros(5,1);
sensorder = ["maxalb","base","minalb"];
misbl = 51;
sensitivitycategories = 0:4;
for bb = 1 : nbblocks
    if preswlk(bb) == false && ~ismember(bb,misbl)
        continue
    elseif ismember(bb,misbl)
        lat1 = block005description(bb,4);
        lon1 = block005description(bb,5);
        aa = areaquad(lat1-10,lon1,lat1,lon1+10,earthellipsoid);
        sensarreas50(1) = sensarreas50(1) + aa;
        continue
    end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    statfilename = strcat(regoutputfiles,"ROstats005_",num2str(bb),".mat");
    load(subfilename,"landmask","pixarea")
    load(radforfname,'GWPmed005')
    
    AOsensmap = ones(size(landmask),'uint8') .* u8noval;    
    AOsensmap(landmask) = 0;
    block_sensarreas50 = zeros(5,1);
    
    if sum(isnan(GWPmed005),'all') == blocksize005^2
        block_sensarreas50(1) = sum(pixarea(landmask),'all','omitnan');
        sensarreas50 = sensarreas50 + block_sensarreas50;
        save(statfilename,"AOsensmap","block_sensarreas50",'-append')
        clear GWPmed005 landmask subfilename AOsensmap statfilename radforfname block_sensarreas50
        continue
    end
    
    load(radforfname,AOarrays{:})
    
    for ss = 1 : numel(sensorder)
        name = strcat("AO",sensorder(ss));
        eval(strcat("map = ",name,";"))
        if ss == 1
            AOsensmap(map < aothreshold) = 4;
        end
        AOsensmap(map >= aothreshold) = numel(sensorder)+1-ss;
    end

    for ss = sensitivitycategories
        kk = AOsensmap == ss;
        block_sensarreas50(ss+1) = sum(pixarea(kk),'all','omitnan');
    end
    sensarreas50 = sensarreas50 + block_sensarreas50;

    
    strcat("Done with creating a sensitivity map in block #",num2str(bb))
    save(statfilename,"AOsensmap","block_sensarreas50",'-append')
    clear(AOarrays{:})
    clear GWPmed005 landmask subfilename AOsensmap statfilename radforfname block_sensarreas50
    clear ss name map ss
end       

save(ReforOppfname,'sensarreas50','-append')