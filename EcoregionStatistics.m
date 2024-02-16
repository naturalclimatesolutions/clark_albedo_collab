% Statistics by Ecoregions


% Calculate Root-to-Shoot Ratios per Ecoregion
% ********************************************
ecor2s = nan(numel(ecoindex)-1,numel(forwsaind)+1);
ecor2sbyindex = nan(numel(ecoindex),numel(forwsaind)+1);
for ee = 1 : numel(ecoindex)
    thisecoid = ecoid(ee);
    [blocks,~] = find(ismember(ecoperblock,thisecoid));
    if ~isempty(blocks) && thisecoid > 0
        r2s = [];forestindexes = [];
        for bb = 1 : numel(blocks)
            thisblock = blocks(bb);
            subfilename = strcat(regoutputfiles,"ROinputs_",num2str(thisblock),".mat");
            load(subfilename,'ecoregion','walkerr2s','modis2010')
            kk = ecoregion == thisecoid;
            vals = walkerr2s(kk);
            r2s = [r2s;vals]; %#ok<AGROW>
            tempind = false(numel(vals),numel(forwsaind));
            for ff = 1 : numel(forwsaind)
                tempind(:,ff) = modis2010(kk) == forwsaind(ff);
            end
            forestindexes = cat(1,forestindexes,tempind);
            clear subfilename kk vals tempind thisblock
        end
        forestindexes = logical(forestindexes);
        for ff = 1 : numel(forwsaind)+1
            if ff == numel(forwsaind)+1
                valmean = mean(r2s,'omitnan');
            else
                kk = forestindexes(:,ff);
                if sum(kk) > 0
                    valmean = mean(r2s(kk),'omitnan');
                else
                    valmean = nan;
                end
            end
            ecor2s(thisecoid,ff) = valmean;
            ecor2sbyindex(ee,ff) = valmean;
        end
        clear valmean r2s tempind forestindexes kk ff
        strcat("Done with calculating root-to-shoot rations in ecoregion #",num2str(thisecoid))
    else
        warning(strcat("Ecoregion #",num2str(thisecoid)," - ",econames(ecoindex(ee)),...
            " is never present in Walker's datamask!"))
    end
    clear thisecoid blocks
end
save(EcoStatfname,'ecor2s','ecor2sbyindex');

% As a first approximation, get the biome means too (even though ecoregions have different sizes)
biomer2s = zeros(numel(biomelist),numel(forwsaind)+1);
for bb = 1 : numel(biomelist)
    ii = biomeid == biomelist(bb);
    biomer2s(bb,:) = mean(ecor2sbyindex(ii,:),1,'omitnan');
end
save(EcoStatfname,'biomer2s','-append');


% Calculate Landcover Areas per Ecoregion or Biome/Region
% *******************************************************
EcoLandcoverAreas = zeros(numel(ecoindex),numel(IGBPBiomes),numel(modisyears));
EcoLandCoverPixelCounts = zeros(numel(ecoindex),numel(IGBPBiomes),numel(modisyears));
BiomesClimatesLancoverAreas = zeros(numel(biomelist),numel(regionnames),numel(kglegnames),...
    numel(IGBPBiomes),numel(modisyears));
BiomesClimatesLandCoverPixelcounts = zeros(numel(biomelist),numel(regionnames),numel(kglegnames),...
    numel(IGBPBiomes),numel(modisyears));

for ee = 1 : numel(ecoindex)
    tic
    thisecoregion = ecoid(ee);
    if thisecoregion == 0, continue; end
    thisbiome = biomeid(ee);
    [blocks,~] = find(ecoperblock==thisecoregion);
    
    if numel(blocks) == 0, continue; end
    
    ecoland = zeros(numel(blocks),numel(IGBPBiomes),numel(modisyears));
    biomeland = zeros(numel(blocks),numel(regionnames),numel(kglegnames),...
        numel(IGBPBiomes),numel(modisyears));
    ecopixs = zeros(numel(blocks),numel(IGBPBiomes),numel(modisyears));
    biomepixs = zeros(numel(blocks),numel(regionnames),numel(kglegnames),...
        numel(IGBPBiomes),numel(modisyears));
    
    for bb = 1 : numel(blocks)
        bid = blocks(bb);
        subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bid),".mat");
        load(subfilename,modisvars{:},'koppengeiger','ecoregion','pixarea','region')
        
        blclim = setdiff(unique(koppengeiger),misskoppengeiger);
        blreg = setdiff(unique(region),regsmis);
        kk = ecoregion == thisecoregion;
        
        for yy = 1 : numel(modisyears)
            eval(strcat("modis = ",modisvars(yy),";"))
            bllc = setdiff(unique(modis(kk)),[waterind,modisnodata]);
            for ll = 1 : numel(bllc)
                lci = bllc(ll);
                ii = modis == bllc(ll); jj = ii+kk==2;
                if sum(isnan(pixarea(jj)),'all') > 0
                    warning(strcat("I have missing pixel areas in block #",num2str(bid),...
                        ", ecoregion #",num2str(thisecoregion)," & landcover: ",IGBPBiomes(lci)))
                end
                ecoland(bb,lci,yy) = sum(pixarea(jj),'all','omitnan');
                ecopixs(bb,lci,yy) = sum(jj,'all');
                for rr = 1 : numel(blreg)
                    oo = region == blreg(rr);
                    for cc = 1 : numel(blclim)
                        mm = koppengeiger == blclim(cc);
                        nn = mm+jj+oo==3;
                        biomeland(bb,blreg(rr),blclim(cc),lci,yy) = ...
                            sum(pixarea(nn),'all','omitnan');
                        biomepixs(bb,blreg(rr),blclim(cc),lci,yy) = sum(nn,'all');
                    end
                end
            end
        end
    end
    EcoLandcoverAreas(ee,:,:) = sum(ecoland,1,'omitnan');
    EcoLandCoverPixelCounts(ee,:,:) = sum(ecopixs,1);
    BiomesClimatesLancoverAreas(thisbiome,:,:,:,:) = BiomesClimatesLancoverAreas(thisbiome,:,:,:,:) + ...
        sum(biomeland,1,'omitnan');
    BiomesClimatesLandCoverPixelcounts(thisbiome,:,:,:,:) = ...
        BiomesClimatesLandCoverPixelcounts(thisbiome,:,:,:,:) + sum(biomepixs,1);

    strcat("Done with statistics in ecoregion #",num2str(thisecoregion),' "',econames(ee),'"')
    clear thisecoregion thisbiome blocks ecoland biomeland bb bid subfilename blclim blreg bllc ...
        kk ll li=ci ii rr oo cc mm nn
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end
save(EcoStatfname,'EcoLandcoverAreas','BiomesClimatesLancoverAreas',...
    'EcoLandCoverPixelCounts','BiomesClimatesLandCoverPixelcounts','-append')


% % Calculate Walker Aboveground Statistics by Ecoregion
% % ****************************************************
% ecowlkagb = nan(numel(ecoindex)-1,numel(despct),2);
% biomewlkagb = nan(numel(biomelist),numel(despct),2);
% ecowlkpix = zeros(numel(ecoindex)-1,2);
% for bio = 1 : numel(biomelist)
%     tic
%     kk = biomeid == bio;
%     ecoidperbio = ecoid(kk);
%     wlkbio = []; reswlkbio = []; clear kk
%     for ee = 1 : numel(ecoidperbio)
%         thisecoid = ecoidperbio(ee);
%         [blocks,~] = find(ismember(ecoperblock,thisecoid));
%         if ~isempty(blocks) && thisecoid > 0
%             wlkecopts = []; wlkecomodispts = [];
%             for bb = 1 : numel(blocks)
%                 thisblock = blocks(bb);
%                 subfilename = strcat(regoutputfiles,"RefOpp005_",num2str(thisblock),".mat");
%                 load(subfilename,'ecoregion','walkeragb','modis2010')
%                 kk = ecoregion == thisecoid;
%                 vals = walkeragb(kk);
%                 vals(vals==0) = nan;
%                 vals = vals(isfinite(vals));
%                 wlkecopts = [wlkecopts;vals]; %#ok<AGROW>
%                 mm = ismember(modis2010,forwsaind);
%                 nn = mm + kk == 2;
%                 resvals = walkeragb(nn);
%                 resvals(resvals==0) = nan;
%                 resvals = resvals(isfinite(resvals));
%                 wlkecomodispts = [wlkecomodispts;resvals]; %#ok<AGROW>
%                 clear thisblock kk vals mm nn resvals
%             end
%             wlkbio = [wlkbio;wlkecopts]; %#ok<AGROW>
%             reswlkbio = [reswlkbio;wlkecomodispts]; %#ok<AGROW>
%             ecowlkagb(thisecoid,:,1) = prctile(wlkecopts,despct);
%             ecowlkagb(thisecoid,:,2) = prctile(wlkecomodispts,despct);
%             ecowlkpix(thisecoid,:) = [numel(wlkecopts),numel(wlkecomodispts)];
%             strcat("Done with calculating Walker AGB statistics in ecoregion #",num2str(thisecoid))
%             clear thisecoid blocks wlkecomodispts wlkecopts bb
%         else
%             warning(strcat("Ecoregion #",num2str(thisecoid)," - ",econames(ecoindex(ee)),...
%                 " is never present in Walker's datamask!"))
%         end
%     end
%     biomewlkagb(bio,:,1) = prctile(wlkbio,despct);
%     biomewlkagb(bio,:,2) = prctile(reswlkbio,despct);
%     clear wlkbio reswlkbio ee kk ecoidperbio
%     strcat("Done with calculating Walker AGB statistics in biome - ",biomenames(bio))
%     eltime = toc;
%     if eltime > 120
%         h = floor(eltime/3600);
%         m = floor((eltime  - h*3600)/60);
%         s = eltime - (h*3600 + m*60);
%         duration([h,m,s])
%     end
% end
% save("I:\Temporary\ReforOpp_add.mat",'ecowlkagb','biomewlkagb','ecowlkpix','ecowlkagb_orig','biomewlkagb_orig');
% 
% 
% % Save other ecoregion characteristics
% % ************************************
% ecoregionsclimates = zeros(numel(ecoindex)-1,numel(kglegnames));
% ecoregionsregions = zeros(numel(ecoindex)-1,numel(regionnames));
% ecoregionslandcover = zeros(numel(ecoindex)-1,numel(IGBPBiomes));
% maxclim = 0; maxreg = 0; maxlc = 0;
% 
% for bb = 1 : nbblocks
%     tic
%     if preswlk(bb) == false, continue; end
%     subfilename = strcat(regoutputfiles,"RefOpp005_",num2str(bb),".mat");
%     load(subfilename,'koppengeiger','ecoregion','region',modisvars{:})
%     
%     dum = ecoperblock(bb,:);
%     ecolist = dum(dum>0);
%     dum = setdiff(unique(ecoregion),[ecomis,0]);
%     if ~isempty(setdiff(dum,ecolist)) || ~isempty(setdiff(ecolist,dum))
%         error(strcat("I have a different set of ecoregions in block #",num2str(bb)))
%     end
%     
%     for ee = 1 : numel(ecolist)
%         thiseco = ecolist(ee);
%         kk = ecoregion == thiseco;
%         
%         climlist = setdiff(unique(koppengeiger(kk)),misskoppengeiger);
%         at = setdiff(ecoregionsclimates(thiseco,:),0);
%         newclimlist = unique([climlist;at']); cl = numel(newclimlist);
%         if cl > maxclim, maxclim = cl; end
%         ecoregionsclimates(thiseco,1:cl) = newclimlist;
%         
%         regionlist = setdiff(unique(region(kk)),regsmis);
%         at = setdiff(ecoregionsregions(thiseco,:),0);
%         newregionlist = unique([regionlist;at']); rg = numel(newregionlist);
%         if rg > maxreg, maxreg = rg; end
%         ecoregionsregions(thiseco,1:rg) = newregionlist;
%         
%         modislist = [];
%         for yy = 1 : numel(modisvars)
%             eval(strcat("modis = ",modisvars(yy),";"))
%             dum = setdiff(unique(modis(kk)),modisnodata);
%             modislist = [modislist;dum]; %#ok<AGROW>
%         end
%         at = setdiff(ecoregionslandcover(thiseco,:),0);
%         newmodislist = unique([modislist;at']); lc = numel(newmodislist);
%         if lc > maxlc, maxlc = lc; end
%         ecoregionslandcover(thiseco,1:lc) = newmodislist;
%         
%     end
%     
%     strcat("done with block #",num2str(bb)," with ",num2str(numel(ecolist))," ecoregions")
% end
% 
% save(ReforOppfname,'ecoregionsclimates','ecoregionsregions','ecoregionslandcover','-append')
% 

