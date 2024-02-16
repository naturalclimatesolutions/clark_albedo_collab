% Import 0.005 maps
%   All Original maps have been projected onto WGS84 and resampled @ 0.005 deg lat-lon, except for
%   Hansen that will be averaged here

% Import Potential Carbon Storage from Walker et al. (2022)
% *******************************
wagb = readgeoraster(WalkerBPotAGBfname);
wbi = georasterinfo(WalkerBPotBGBfname);
wbgb = readgeoraster(WalkerBPotBGBfname);
if wbi.RasterSize(1) ~= walkerinfo.RasterSize(1) || ...
        wbi.RasterSize(2) ~= walkerinfo.RasterSize(2)
    wbgb = snaparray(walkerinfo,wbi,wbgb);
end
wlkgeodat = walkerinfo.RasterReference.GeographicCRS;
wlkmis = walkerinfo.MissingDataIndicator;

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    y = ceil(bb/nnblx);
    x = bb - (y-1)*nnblx;
    smregi = firstiW005(y) : firstiW005(y) + blocksize005 - 1;
    smregj = firstjW005(x) : firstjW005(x) + blocksize005 - 1;
    
    northlat = ncs(y); southlat = ncs(y+1);
    westlon = wcs(x); eastlon = wcs(x+1);
    blockinfo = georefcells([southlat northlat],[westlon eastlon],...
        [blocksize005 blocksize005],'ColumnsStartFrom','north');
    blockinfo.GeographicCRS = wlkgeodat;

    walkeragb = double(wagb(smregi,smregj));
    walkeragb(walkeragb<0) = nan;
    walkeragb(walkeragb==wlkmis) = nan;
    if isempty(find(isfinite(walkeragb),1))
        preswlk(bb) = false; %#ok<SAGROW>
        strcat("No data of Walker potential carbon in block #",num2str(bb))
    else
        walkerbgb = double(wbgb(smregi,smregj));
        walkerbgb(walkerbgb<0) = nan;
        walkerr2s = walkerbgb ./ walkeragb;
        walkerr2s(isinf(walkerr2s)) = nan;
        walkertotbio = walkeragb + walkerbgb;
        gapfillwalker = false;
        if sum(isfinite(walkeragb),'all') ~= sum(isfinite(walkerbgb),'all')
            gapfillwalker = true;
        end
        walkerdatamask = walkertotbio >= 0;
        if sum(walkerdatamask,'all') == 0
            preswlk(bb) = false;%#ok<SAGROW>
            strcat("No data of Walker potential carbon in block #",num2str(bb))
        else
            save(subfilename,'x','y','smregi','smregj','northlat','southlat','westlon','eastlon',...
                'blockinfo','walkeragb','walkerbgb','walkerr2s','walkertotbio',...
                'walkerdatamask','gapfillwalker')
            strcat("Done with importing Walker potential carbon into block #",num2str(bb))
        end

    end
    
    clear x y walkeragb walkerbgb walkerr2s smregi smregj norhtlat southlat eastlon westlon ...
        blockinfo subfilename walkertotbio walkerdatamask gapfillwalker 

end
save(sparameterfile,'wlkgeodat','preswlk','wlkmis','-append')
clear wagb wbgb


% Import MODIS land cover (MCD12Q1 2001 and 2010)
% ***********************
modisperblock = zeros(10,8,nbblocks);
modisfilesperblock = strings(nbblocks,10);
modisvars = strings(numel(modisyears),1);
for yy = 1 : numel(modisyears)
    modisvars(yy) = strcat("modis",modisyears(yy));
end
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'northlat','southlat','westlon','eastlon','blockinfo','walkerdatamask')
    
    ll = logical(ismember(table2array(modiscoord(:,strcmp(modisccvars,"lat_max"))),northlat) + ...
        ismember(table2array(modiscoord(:,strcmp(modisccvars,"lat_min"))),southlat));
    kk = (table2array(modiscoord(:,strcmp(modisccvars,"lon_min"))) <= westlon) + ...
        (table2array(modiscoord(:,strcmp(modisccvars,"lon_max"))) >= westlon) == 2;
    nn = (table2array(modiscoord(:,strcmp(modisccvars,"lon_max"))) >= eastlon) + ...
        (table2array(modiscoord(:,strcmp(modisccvars,"lon_min"))) <= eastlon) == 2;
    
    oo = logical(kk+nn);
    mod = modiscoord(ll+oo==2,:); clear ll kk nn
    [nf,~] = size(mod);
    modisfiles = strings(nf,2);
    for kk = 1 : nf
        hp = num2str(table2array(mod(kk,strcmp(modisccvars,"ih"))));
        if length(hp) == 1, hp = strcat('0',hp); end
        vp = num2str(table2array(mod(kk,strcmp(modisccvars,"iv"))));
        if length(vp) == 1, vp = strcat('0',vp); end
        % check if the modis file actually exists
        modtile = strcat("h",hp,"v",vp);
        if sum(contains(modis2001filelist005,modtile)) == 1 && ...
                sum(contains(modis2010filelist005,modtile)) == 1
            modisfiles(kk,1) = strcat("MLC01.",modtile,".tif");
            modisfiles(kk,2) = strcat("MLC10.",modtile,".tif");
        elseif sum(contains(modis2001filelist005,modtile)) == 1 || ...
                sum(contains(modis2010filelist005,modtile)) == 1
            error(strcat("There is discrepancy between 2001 and 2010 for ",modtile))
        elseif sum(contains(modisoriginalfilelist,modtile)) == 1
            error(strcat('MODIS file with tile "',modtile,'" is missing for 2001 & 2010 in folder "',...
                modis005folder,'"'))
        end
    end
    dd = modisfiles(:,1) == ""; vv = logical(1-dd);
    modisfiles = modisfiles(vv,:); mod = mod(vv,:); nf = nf - sum(dd);
    modisfilesperblock(bb,1:nf) = modisfiles(:,1);

    % check "real" longitude limits as they might be different from original due to projection
    if nf > 0        
        mfl = true(nf,1);
        lon = zeros(nf,2);
        lat = zeros(nf,2);
        for ff = 1 : nf
            mfname = strcat(modis005folder,"2001\",modisfiles(ff,1));
            modisinfo = georasterinfo(mfname);
            lon(ff,:) = modisinfo.RasterReference.LongitudeLimits;
            lat(ff,:) = modisinfo.RasterReference.LatitudeLimits;
            if lon(ff,2) <= westlon || lon(ff,1) >= eastlon
                mfl(ff) = false;
            end
            clear mfname modisinfo dum
        end
        modisfiles = modisfiles(mfl,:); nf = sum(mfl);
        if exist('mod','var') == 1
            mod = mod(mfl,:); lat = lat(mfl,:); lon = lon(mfl,:);
            modisperblock(1,1:4,bb) = [westlon,eastlon,southlat,northlat];
            modisperblock(2:nf+1,1:4,bb) = table2array(mod(:,3:6));
            modisperblock(2:nf+1,5:6,bb) = lon;
            modisperblock(2:nf+1,7:8,bb) = lat;
        end
        for yy = 1 : numel(modisyears)
            year = modisyears(yy);
            mfname = strcat(modis005folder,year,"\",modisfiles(1,yy));
            modisinfo = georasterinfo(mfname);
            modisnodata = modisinfo.MissingDataIndicator;
            modisval = readgeoraster(mfname);
            if modisinfo.RasterSize(1) ~= blockinfo.RasterSize(1) || ...
                    modisinfo.RasterSize(2) ~= blockinfo.RasterSize(2)
                modisval = snaparray(blockinfo,modisinfo,modisval);
            end
            if nf > 1
                for mm = 2 : nf
                    mfname = strcat(modis005folder,year,"\",modisfiles(mm,yy));
                    modisinfo = georasterinfo(mfname);
                    modisval2 = readgeoraster(mfname);
                    if modisinfo.RasterSize(1) ~= blockinfo.RasterSize(1) || ...
                            modisinfo.RasterSize(2) ~= blockinfo.RasterSize(2)
                        modisval2 = snaparray(blockinfo,modisinfo,modisval2);
                    end
                    modisval(modisval==modisnodata) = modisval2(modisval==modisnodata);
                    clear mfname modisinfo modisval2
                end
            end

            % Because of mismatch of projections, I have some pixels missing between two adjacent
            % modis tiles, therefore, I will have to gapfill here, but on a small search radius
            fakemask = true(size(modisval));
            fakemask(modisval==17) = false;
            modisval = blockgapfill(modisval,fakemask,blocksize005,modisnodata,5);
            if sum(modisval==modisnodata,'all') > 0
                if sum(modisval(walkerdatamask) == modisnodata,'all') > 0
                    modisval = blockgapfill(modisval,walkerdatamask,blocksize005,modisnodata,50);
                    npts = sum(modisval(walkerdatamask) == modisnodata,'all');
                    if npts > 0
                        warning(strcat("There are ",num2str(npts)," missing MODIS ",...
                            num2str(modisyears(yy))," values in block #",num2str(bb)))
                    end
                end
            end
            mask = ismember(modisval,landind);
            if yy == 1
                modislandmask = mask;
            else
                modislandmask = logical(modislandmask + mask);
            end
            eval(strcat(modisvars(yy)," = modisval;"))
        end
        save(subfilename,modisvars{:},'modislandmask','-append')
    else
        if northlat >= 0, lat = strcat(num2str(northlat),"N"); else, lat = strcat(num2str(abs(northlat)),"S"); end
        if westlon >= 0, lon = strcat(num2str(westlon),"E"); else, lon = strcat(num2str(abs(westlon)),"W"); end
        warning(strcat("there is no MODIS landcover values for block #",num2str(bb),...
            " [",lat,",",lon,"]"))
    end
    
    strcat("done with block #",num2str(bb)," - with ",num2str(nf)," MODIS files")
    clear subfilename northlat southlat westlon eastlon bioinfo ll kk nn oo mod nf modisfiles kk ...
        hp vp mfname modisinfo modisval mm modisval2 dd
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end
save(sparameterfile,'modisperblock','modisfilesperblock','modisvars','modisnodata','-append')
clear(modisvars{:})


% Consolidate land masks
% **********************
% (Land mask is determined by modis, unless some Walker pixels need "land")
for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'modislandmask','walkerdatamask')
    if exist('modislandmask','var') == 0
        landmask = walkerdatamask;
    else
        landmask = modislandmask;
    end
    landmask(walkerdatamask) = true;
    save(subfilename,'landmask','-append')
end


% Import Ecoregions from Dinerstein et al. (2017)
% *****************
ecoinfo = georasterinfo(ecoregion005fname);
ecoval = readgeoraster(ecoregion005fname);
if ecoinfo.RasterReference.RasterSize(1) ~= walkerinfo.RasterReference.RasterSize(1) || ...
        ecoinfo.RasterReference.RasterSize(2) ~= walkerinfo.RasterReference.RasterSize(2)
    ecoval = snaparray(walkerinfo,ecoinfo,ecoval);
end
ecomis = ecoinfo.MissingDataIndicator;
ecoperblock = nan(nbblocks,100); nbmaxeco = 0;
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj','landmask')
    ecoregion = ecoval(smregi,smregj);
    ecoregion = blockgapfill(ecoregion,landmask,blocksize005,ecomis,50);  % gap-fill
    blecolst = setdiff(unique(ecoregion),ecomis);
    ecoperblock(bb,1:numel(blecolst)) = blecolst;
    if numel(blecolst) > nbmaxeco; nbmaxeco = numel(blecolst); end
    save(subfilename,'ecoregion','-append')
    clear subfilename ecoregion blecolst
    strcat("Done with importing Dinerstein Ecoregion in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
end
ecoperblock = ecoperblock(:,1:nbmaxeco);
save(sparameterfile,'ecomis','ecoperblock','-append')
clear ecoval


% Import 25-regions
% *****************
regsinfo = georasterinfo(reg005fname);
regsval = readgeoraster(reg005fname);
if regsinfo.RasterSize(1) ~= walkerinfo.RasterSize(1) || ...
        regsinfo.RasterSize(2) ~= walkerinfo.RasterSize(2)
    regsval = snaparray(walkerinfo,regsinfo,regsval);
end
regsmis = setdiff(unique(regsval),1:25);

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj','landmask');
    
    region = regsval(smregi,smregj);
    region = blockgapfill(region,landmask,blocksize005,regsmis,50);  % gap-fill

    save(subfilename,'region','-append');
    
    clear smregi smregj region subfilename walkerdatamask
    
end
save(sparameterfile,'regsinfo','regsmis','-append')
clear regsval


% Import Restoration categories from Walker et al. (2022)
% *****************************
wlkocatfname = strcat(walkerfolder,"Walker_OCat_WGS84.tif");
wlkocat = readgeoraster(wlkocatfname);
walkercategoriesinfo = georasterinfo(wlkocatfname);
if walkercategoriesinfo.RasterSize(1) ~= walkerinfo.RasterSize(1) || ...
        walkercategoriesinfo.RasterSize(2) ~= walkerinfo.RasterSize(2)
    wlkocat = snaparray(walkerinfo,walkercategoriesinfo,wlkocat);
end
wlkocmis = walkercategoriesinfo.MissingDataIndicator;

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj');
    
    walkeroppcat = wlkocat(smregi,smregj);
    save(subfilename,'walkeroppcat','-append');
    
    clear smregi smregj walkeroppcat subfilename
    
end
save(sparameterfile,'walkercategoriesinfo','wlkocmis','-append')
clear wlkocat


% Import Tree restoration potential from Bastin et al. (2019)
% *********************************
bastininfo = georasterinfo(Bastinfname);
btval = readgeoraster(Bastinfname);
bastinmiss = nan;
rs = bastininfo.RasterSize;
g = walkerinfo.RasterReference.GeographicCRS;
% somehow Matlab cannot read the info of the file properly, and because the lat/lon limits are
% outside of the normal range, I will just delete the "offending" cells manually!
corbtval = btval(3:rs(1),3:rs(2));
bastininfo = georasterref('RasterSize',size(btval),'RasterInterpretation','cells',...
      'ColumnsStartFrom','north','LatitudeLimits', [-60.000075 90.009925-2*latlon005],...
       'LongitudeLimits', [-180.009850+2*latlon005 180.000150]);
bastininfo.GeographicCRS = g;

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj');
    
    bastin = corbtval(smregi,smregj);
    save(subfilename,'bastin','-append');
    
    clear smregi smregj bastin subfilename
end
save(sparameterfile,'bastininfo','bastinmiss','-append')
clear bastininfo btval


% Import Total Potential Tree Cover from Bastin et al. (2019)
% *********************************
bastinTPinfo = georasterinfo(BastinTPfname);
btval = readgeoraster(BastinTPfname);
bastinTPmiss = bastinTPinfo.MissingDataIndicator;
if bastinTPinfo.RasterSize(1) ~= walkerinfo.RasterSize(1) || ...
        bastinTPinfo.RasterSize(2) ~= walkerinfo.RasterSize(2)
    btval = snaparray(walkerinfo,bastinTPinfo,btval);
end

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj');
    
    bastinTP = btval(smregi,smregj);
    save(subfilename,'bastinTP','-append');
    
    clear smregi smregj bastinTP subfilename
end
save(sparameterfile,'bastinTPinfo','bastinTPmiss','-append')
clear bastinTPinfo btval


% Import Natural Climate Solutions from Griscom et al. (2017)
% ********************************
griscominfo = georasterinfo(Griscomfname);
gcval = readgeoraster(Griscomfname);
if griscominfo.RasterReference.RasterSize(1) ~= walkerinfo.RasterReference.RasterSize(1) || ...
        griscominfo.RasterReference.RasterSize(2) ~= walkerinfo.RasterReference.RasterSize(2)
    gcval = snaparray(walkerinfo,griscominfo,gcval);
end
missgriscom = griscominfo.MissingDataIndicator;

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj');
    
    griscom = gcval(smregi,smregj);
    save(subfilename,'griscom','-append');
    
    clear smregi smregj griscom subfilename
end
save(sparameterfile,'griscominfo','missgriscom','-append')
clear gcval griscominfo


% Import Koppen-Geiger climate map
% ********************************
koppengeigerinfo = georasterinfo(KGClimatefname);
kgval = readgeoraster(KGClimatefname);
if koppengeigerinfo.RasterReference.RasterSize(1) ~= walkerinfo.RasterReference.RasterSize(1) || ...
        koppengeigerinfo.RasterReference.RasterSize(2) ~= walkerinfo.RasterReference.RasterSize(2)
    kgval = snaparray(walkerinfo,koppengeigerinfo,kgval);
end
misskoppengeiger = setdiff(unique(kgval),climateid);

for bb = 1 : nbblocks
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'smregi','smregj');
        
    koppengeiger = kgval(smregi,smregj);
    koppengeiger = blockgapfill(koppengeiger,modislandmask,blocksize005,misskoppengeiger,50);  % gap-fill
    save(subfilename,'koppengeiger','-append');
    
    clear x y smregi smregj koppengeiger subfilename
end
save(sparameterfile,'koppengeigerinfo','misskoppengeiger','-append')
clear koppengeigerinfo kgval


% Import Hansen Tree Cover (from orginal data Hansen et. al (2013)
% ************************
filepref = strcat(origdatafolder,"HansenTreeCover\TreeCover\Hansen_GFC-2019-v1.7_treecover2000_");
missfiles = [];
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    subfilename = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(subfilename,'northlat','westlon')
    if northlat == 0
        dum1 = "00N";
    elseif northlat > 0
        dum1 = strcat(num2str(northlat),'N');    
    else
        dum1 = strcat(num2str(abs(northlat)),'S');
    end
    if westlon == 0
        dum2 = "000E";
    elseif westlon < 0
        lonv = abs(westlon);
        if lonv >= 100, dum2 = strcat(num2str(lonv),'W'); ...
        else, dum2 = strcat('0',num2str(lonv),'W'); end
    else
        if westlon >= 100, dum2 = strcat(num2str(westlon),'E'); ...
        else, dum2 = strcat('0',num2str(westlon),'E'); end
    end
    filesuff = strcat(dum1,"_",dum2,".tif");
    fname = strcat(filepref,filesuff);
    if exist(fname,'file') == 0
        missfiles = cat(1,missfiles,filesuff);
        strcat("Missing Hansen Tree Cover data for block #",num2str(bb))
    else
        treecover = zeros(blocksize005,blocksize005);
        orig = readgeoraster(fname);
        hts = size(orig);
        tilesize = hts(1) ./ blocksize005;        
        for i = 1 : blocksize005
            ii = (i-1)*tilesize+1 : i*tilesize;
            for j = 1 : blocksize005
                jj = (j-1)*tilesize+1 : j*tilesize;
                treecover(i,j) = mean(orig(ii,jj),'all');
            end
        end
        save(subfilename,'treecover','-append')
    end

    strcat("Done with importing Hansen Tree Cover in block #",num2str(bb))
    clear subfilename northlat westlon dum1 dum2 lonv filesuff fname treecover orig hts ...
        i ii j jj tilesize
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
end
   
    
