% NCS - Global Albedo - TNC-Clark Collaboration - Combine Albedo & Carbon

% Gap-fill Walker BGB (if needed)
% *******************
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'gapfillwalker')
    if gapfillwalker == true
        load(subfilename,'walkeragb','walkerbgb','walkerr2s','modis2010',...
            'ecoregion','walkerdatamask')
        walkerbgb(walkeragb==0) = 0;    % have zeros instead of nan later        
        aa = isfinite(walkeragb);
        cc = isfinite(walkerbgb);
        ff = (aa + cc) == 1;
        [i,j] = find(ff);
        for k = 1 : numel(i)
            ii = i(k);
            jj = j(k);
            pteco = ecoregion(ii,jj); 
            if pteco == ecomis
                warning(strcat("No ecoregion ID here! (in block #",num2str(bb),")"))
            else
                ee = ecoindex == pteco;
                ptforest = modis2010(ii,jj); ff = forwsaind == ptforest;
                if sum(ff) == 1
                    r2s = ecor2s(ee,ff,despct == 50);
                else
                    dum = squeeze(ecor2s(ee,:,despct == 50)); dum(dum==0) = nan;
                    r2s = mean(dum,'omitnan');
                end
                if isnan(walkerbgb)
                    walkerbgb(ii,jj) = walkeragb(ii,jj) .* r2s;
                else
                    walkeragb(ii,jj) = walkerbgb(ii,jj) ./ r2s;
                end
                walkerr2s(ii,jj) = r2s;
                walkerdatamask(ii,jj) = true;
            end
        end
        walkertotbio = walkeragb + walkerbgb;
        save(subfilename,'walkeragb','walkerbgb','walkertotbio','walkerr2s',...
            'walkerdatamask','-append')
        strcat("Done gap-filling walker in block #",num2str(bb))
    end
end


% "Gap-fill" pixelareas
% *********************
Wsc = WlkClassID(contains(WlkClassDef,"R/"));
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    load(subfilename,"selectedbastin","walkeroppcat","griscom","pixarea",...
        'northlat','southlat','eastlon','westlon')
    load(radforfname,"RFmed005")
    
    latser = northlat : - latlon005 : southlat;
    lonser = westlon : latlon005 : eastlon;

    mask1 = isfinite(RFmed005);
    mask2 = ismember(walkeroppcat,Wsc);
    mask3 = griscom == 1;
    datmask = logical(mask1 + mask2 + mask3 + selectedbastin);
    dum = datmask - isfinite(pixarea); misarea = dum == 1;

    [a,b] = find(misarea);
    for pt = 1 : numel(a)
        ii = a(pt);
        jj = b(pt);
        lat1 = latser(ii); lat2 = latser(ii+1);
        lon1 = lonser(jj); lon2 = lonser(jj+1);
        pixarea(ii,jj) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
    end

    if numel(a) > 0
        save(subfilename,'pixarea','-append')
        strcat("Done gap-filling pixel area in block #",num2str(bb))
    end
    clear("RFmed005","selectedbastin","walkeroppcat","griscom","pixarea",...
        'northlat','southlat','eastlon','westlon')
    clear mask1 mask2 mask3 damask dum misarea a b pt ii jj lat1 lon1 lat2 lon2
end
save(sparameterfile,"Wsc",'-append')


% Create Albedo Offset and Net Climate Impact maps (no stats)
% ************************************************
albedostats = ["med","min","max"];
bdnames = ["base","minalb","maxalb"];   % these will be minimum and maximum effect, which will differ from RFmin and RFmax
AOarrays = strcat("newAO",bdnames);
NCIarrays = strcat("newNCI",bdnames);
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    load(subfilename,"walkertotbio")
    load(radforfname,GWPmaps{:})
    
    if sum(isnan(GWPmax005),'all') ~= blocksize005^2
        wlk = walkertotbio * 44/12;
        
        for ss = 1 : numel(albedostats)
            eval(strcat("alb = ",GWPmaps(ss),";"));
            AO = (-alb ./ wlk) .* 100;
%             posrf = alb>0;
%             kinf = isinf(AO);       % isinf does not differenciate -inf and inf ...
%             kneginf = kinf; kneginf(logical(1-posrf)) = false;
%             kinf(posrf) = false;
%             AO(kinf) = offsetcat(numel(offsetcat));
%             AO(kneginf) = offsetcat(1);

            NCI = alb + wlk;

            eval(strcat(AOarrays(ss)," = AO;"))
            eval(strcat(NCIarrays(ss)," = NCI;"))
            clear alb AO posrf kinf kneginf NCI
        end
        WalkertotCO2 = wlk;
        WlkTCO2RF = wlk;
        WlkTCO2RF(isnan(GWPmed005)) = nan;
        save(radforfname,AOarrays{:},NCIarrays{:},"WalkertotCO2","WlkTCO2RF",'-append')
        clear(AOarrays{:},NCIarrays{:})
        strcat("Done with calculating low/high bounds of AO and NCI in block #",num2str(bb))
    else
        strcat("No valid values of albedo in block #",num2str(bb))
    end
    clear("walkertotbio",GWPmaps{:})
end

save(sparameterfile,"bdnames","AOarrays","NCIarrays",'-append')








