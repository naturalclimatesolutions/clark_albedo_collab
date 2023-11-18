% Truncate Walker, calculate new AO and NCI
% *****************************************
trint = permute(despct(despct>75),[2,1]);
nmaps = numel(trint);
newagb = strcat("walker",num2str(trint));
walkeresatrunc = strcat("walkertot",num2str(trint));
AOesatrmaps = strcat("AOwlk",num2str(trint));
NCIesatrmaps = strcat("NCIwlk",num2str(trint));
reducedesa = esaperecoid(:,:,ismember(despct,trint));
aolowcap = offsetcat(1);
aohighcap = offsetcat(numel(offsetcat));
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    load(subfilename,"walkerdatamask")
    load(radforfname,"GWPmed005","WalkertotCO2")
    if sum(isnan(GWPmed005),'all') == blocksize005^2
        aoflag = false;
    else
        aoflag = true;
        alb = GWPmed005;
    end
    clear GWPmed005
    if presbio(bb) == false
        map = nan(blocksize005,blocksize005,'single');
        map(walkerdatamask) = 0;
        for mm = 1 : nmaps
            eval(strcat(walkeresatrunc(mm)," = map;"))
        end
        regpix = sum(isfinite(WalkertotCO2),'all');
        truncpix = sum(isfinite(walkertot85),'all');
        if regpix ~= truncpix
            error(strcat("My pixel count is not the same in block #"),num2str(bb))
        end

        save(radforfname,walkeresatrunc{:},'-append')
        if aoflag == true
            ao = nan(blocksize005,blocksize005,'single');
            ao(alb>0) = aolowcap;
            ao(alb<0) = aohighcap;
            ao(walkerdatamask==0) = nan;
            nci = alb;
            nci(walkerdatamask==0) = nan;
            for mm = 1 : nmaps
                eval(strcat(AOesatrmaps(mm)," = ao;"))
                eval(strcat(NCIesatrmaps(mm)," = nci;"))
            end
            save(radforfname,AOesatrmaps{:},NCIesatrmaps{:},'-append')
            strcat("Walker-truncated maps (no ESA but AO/NCI) in block #",num2str(bb))
        else
            strcat("Walker-truncated maps (no ESA/AO/NCI) in block #",num2str(bb))
        end
        clear(walkeresatrunc{:},"walkerdatamask")
    else
        load(subfilename,'walkeragb','walkerr2s','walkerbgb',"ecoregion")
        load(strcat(regoutputfiles,"MostLikely005_",num2str(bb),".mat"),"MLFC");
        
        r2s = walkerr2s;
        r2s(walkeragb==0) = 0;
        % I corrected the gap-fill in CombineAlbedoCarbon now, but did not
        % rerun it, so this would be obsolete in the future
        prestot = isfinite(WalkertotCO2);
        if sum(isnan(r2s(prestot))) > 0
            pxs = find(prestot + isnan(r2s) == 2);
            for k = 1 : numel(pxs)
                r2s(pxs(k)) = walkerbgb(pxs(k)) ./ walkeragb(pxs(k));
            end
            if sum(isnan(r2s(prestot))) > 0 
                error(strcat("Not sure how I still have nan values in r2s where tot bio is finite",...
                    " n block # ",num2str(bb)))
            end
        end
        agbw = single(walkeragb);
        thiseco = setdiff(unique(ecoregion),[0,u16noval]);
        truncwalker = repmat(agbw,[1,1,nmaps]);
        if aoflag == true
            AOtrunc = truncwalker;
            NCItrunc = truncwalker;
        end
        
        parfor mm = 1 : nmaps
            % create truncated Walker
            agb = truncwalker(:,:,mm);
            esastats = reducedesa(:,:,mm);
            for ee = 1 : numel(thiseco)
                eco = thiseco(ee);
                ecopix = ecoregion == eco;
                fteco = unique(MLFC(ecopix)); %#ok<PFBNS>
                for ff = 1 : numel(fteco)
                    ft = fteco(ff);
                    ffpix = MLFC == ft;
                    if ismember(ft,[0,u8noval]), fti = numel(forwsaind) + 1; else, fti = find(forwsaind==ft); end
                    val = esastats(eco,fti);
                    if isnan(val), val = esastats(eco,numel(forwsaind)+1); end
                    if ~isnan(val)
                        vpix = walkeragb > val;
                        kk = ecopix + ffpix + vpix == 3;
                        agb(kk) = val;
                    else
                        kk = ecopix + ffpix == 2;
                        agb(kk) = 0;
                    end
                end
            end
            tot = agb .* (1 + r2s) .* 44/12;
            if ~isempty(find(tot-WalkertotCO2>10^-4,1))
                pts = find(tot-WalkertotCO2>10^-4);
                error(strcat("I have ",num2str(numel(pts))," pixels where ",num2str(trint(mm)),...
                    "% truncated values are above original ones in block #",num2str(bb)))
            end
            pts = find(tot-WalkertotCO2>0);
            tot(pts) = WalkertotCO2(pts);
            if sum(isnan(tot),'all') ~= sum(isnan(WalkertotCO2),'all')
                if abs(sum(isnan(tot),'all') - sum(isnan(WalkertotCO2),'all')) / ...
                       sum(isnan(WalkertotCO2),'all') > 0.001
                error(strcat("check nans in block #",num2str(bb)))
                end
            end

            % Calculate AO and NCI where available
            if aoflag == true
                ao = (-alb ./ tot) .* 100;
                ao(ao<offsetcat(1)) = offsetcat(1);
                ao(ao>offsetcat(numel(offsetcat))) = offsetcat(numel(offsetcat));
%                 posrf = alb>0;
%                 kinf = isinf(ao);       % isinf does not differenciate -inf and inf ...
%                 kneginf = kinf; kneginf(logical(1-posrf)) = false;
%                 kinf(posrf) = false;
%                 ao(kinf) = infvals(2);    %#ok<PFBNS>
%                 ao(kneginf) = infvals(1);
                nci = alb + tot;
                AOtrunc(:,:,mm) = ao;
                NCItrunc(:,:,mm) = nci;
            end
            truncwalker(:,:,mm) = tot;
        end
        
        for mm = 1 : nmaps
            eval(strcat(walkeresatrunc(mm)," = truncwalker(:,:,mm);"))
            if aoflag == true
                eval(strcat(AOesatrmaps(mm)," = AOtrunc(:,:,mm);"))
                eval(strcat(NCIesatrmaps(mm)," = NCItrunc(:,:,mm);"))
            end
        end

        regpix = sum(isfinite(WalkertotCO2),'all');
        truncpix = sum(isfinite(walkertot85),'all');
        if regpix ~= truncpix
            error(strcat("My pixel count is not the same in block #",num2str(bb)))
        end
        if aoflag == true
            save(radforfname,walkeresatrunc{:},AOesatrmaps{:},NCIesatrmaps{:},'-append')
        else
            save(radforfname,walkeresatrunc{:},'-append')
        end
        clear(walkeresatrunc{:},AOesatrmaps{:},NCIesatrmaps{:},'walkeragb','walkerr2s',"ecoregion","MLFC")
        clear truncwalker AOtrunc NCItrunc mm aoflag alb origwalker thiseco
        
        strcat("Done with creating Walker-truncated maps and AO/NCI in block #",num2str(bb))
        eltime = toc;
        if eltime > 120
            h = floor(eltime/3600);
            m = floor((eltime  - h*3600)/60);
            s = eltime - (h*3600 + m*60);
            duration([h,m,s])
        end
    end
end

save(sparameterfile,"walkeresatrunc","AOesatrmaps","NCIesatrmaps",'-append')