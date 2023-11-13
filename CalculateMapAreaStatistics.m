% Calculate area statistics of Carbon, Albedo, NCI and AO (for stats bar in map figure)
% *******************************************************
% In this version I will limit to the maps we really need, which are:
%   - median Radiative Forcing globally
%   - Walker AGB+BGB in CO2e globally (limited to pixels with RF)
%   - Albedo Offsets globally (median RF)
%   - Net Climate impact globally and by reforestation opportunities (median RF)
% Note that I will 

cvals = -2000:2000; ncvals = numel(cvals);
aovals = -10000:10000; naovals = numel(aovals);
redarr = ismember(aovals,cvals);
ncilist = strcat(allreforestationopp,"NCI");
basename = cat(1,["AlbedoGWP";"AlbedoOffset";"CarbonOnlyRFlimited"],permute(ncilist,[2,1]));
areacountarraylist = strcat(basename,"_AreaCount");
varlist = ["GWPmed005","AObase","WlkTCO2RF","NCIbase","NCIbase","NCIbase","NCIbase"];
blockarrays = strcat("Block_",areacountarraylist);
pfsize = [blocksize005,blocksize005,numel(areacountarraylist)];

for aa = 1 : numel(areacountarraylist)
    aname = areacountarraylist(aa);
    nvals = ncvals;
    if contains(aname,"AlbedoOffset"), nvals = naovals; end
    eval(strcat(aname," = zeros(nvals,1);"))
end


% Sum areas for each value
% ************************
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    load(subfilename,"landmask","pixarea")
    load(radforfname,"GWPmed005","WlkTCO2RF")
    statfilename = strcat(regoutputfiles,"ROstats005_",num2str(bb),".mat");

    if sum(isnan(GWPmed005),'all') == blocksize005^2, continue; end % To compare maps, limit Walker to RF pixels
    load(radforfname,"AObase","NCIbase") %#ok<NASGU>
    load(subfilename,"selectedbastin","walkeroppcat","griscom","combinedopp")
    masks = false(pfsize);
    data = nan(pfsize,'single');
    intarray = zeros(naovals,numel(areacountarraylist));

    for aa = 1 : numel(areacountarraylist)
        bname = basename(aa);
        eval(strcat("data(:,:,aa) = round(",varlist(aa),");"))

        switch bname
            case "WalkerNCI"
                masks(:,:,aa) = ismember(walkeroppcat,Wsc);
            case "GriscomNCI"
                masks(:,:,aa) = griscom == 1;
            case "BastinNCI"
                masks(:,:,aa) = selectedbastin;
            case "CombinedOppNCI"
                masks(:,:,aa) = combinedopp;
            otherwise
                masks(:,:,aa) = landmask;
        end
        clear aname
    end

    parfor aa = 1 : numel(areacountarraylist)
        ardata = data(:,:,aa);
        ararray = intarray(:,aa);
        armask = masks(:,:,aa);
        ardata(armask==0) = nan;
        ardata(ardata<aovals(1)) = aovals(1); %#ok<PFBNS> 
        ardata(ardata>aovals(numel(aovals))) = aovals(numel(aovals));
        arpixarea = pixarea;
        datalist = unique(ardata(isfinite(ardata)));
        for cc = 1 : numel(datalist)
            cid = aovals == datalist(cc);
            pix = ardata == datalist(cc);
            ararray(cid) = sum(arpixarea(pix),'all','omitnan');
        end
        intarray(:,aa) = ararray;
    end

    for aa = 1 : numel(areacountarraylist)
        aname = blockarrays(aa);
        array = intarray(:,aa);
        if ~contains(aname,"AlbedoOffset"), array = array(redarr); end
        eval(strcat(aname," = array;"))
        eval(strcat(areacountarraylist(aa)," = ",areacountarraylist(aa)," + array;"))
        clear aname array
    end
    clear aa masks data intarray

    save(statfilename,blockarrays{:})
    clear(blockarrays{:},"WalkertotCO2","AObase","NCIbase","GWPmed005","selectedbastin",...
        "walkeroppcat","griscom","landmask","pixarea")

    strcat("Done with adding areas for map stats in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
end


% Calculate statistics
% ********************
statsarraylist = strcat(areacountarraylist,"Stats");
for aa = 1 : numel(areacountarraylist)
    eval(strcat("data = ",areacountarraylist(aa),";"))
    vals = cvals;
    if contains(areacountarraylist(aa),"AlbedoOffset"), vals = aovals; end
    statsarray = zeros(numel(mappct),1);
    totarea = sum(data);
    cumdata = cumsum(data);
    areastval = totarea .* mappct / 100;
    for ss = 1 : numel(mappct)
        ii = find(cumdata>=areastval(ss),1);
        statsarray(ss) = vals(ii); clear ii
    end
    eval(strcat(statsarraylist(aa)," = statsarray;"))
    clear data vals statsarray totarea cumdata areastval ss
end

save(ReforOppfname, statsarraylist{:},areacountarraylist{:})
save(sparameterfile,"areacountarraylist","statsarraylist",'-append')