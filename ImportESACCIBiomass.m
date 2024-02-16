% Import ESA-CCI and calculate ecoregion statistics at the original resolution of 100m
% ************************************************************************************

% Get the high resolution biomass values from ESA-CCI (and save in mat-files)
% ***************************************************
for bb = 1 : nbblocks
    subfilename = strcat(subdata10,"BioSubData0008_",num2str(bb),".mat");
    y = ceil(bb/nblx);
    x = bb - (y-1)*nblx;
    northlat = ncs(y); southlat = ncs(y+1);
    westlon = wcs(x); eastlon = wcs(x+1);
    if northlat == 0
        dum1 = "N00";
    elseif northlat > 0
        dum1 = strcat('N',num2str(northlat));    
    else
        dum1 = strcat('S',num2str(abs(northlat)));
    end
    if westlon == 0
        dum2 = "E000";
    elseif westlon < 0
        lonv = abs(westlon);
        if lonv >= 100, dum2 = strcat('W',num2str(lonv)); else, dum2 = strcat('W0',num2str(lonv)); end
    else
        if westlon >= 100, dum2 = strcat('E',num2str(westlon)); ...
        else, dum2 = strcat('E0',num2str(westlon)); end
    end
    coordname = strcat(dum1,dum2);
    bioline = contains(biofilelist,coordname);
    if sum(bioline) == 1
        bioname = strcat(biomassfolder,biofilelist(bioline));
        bioinfo = georasterinfo(bioname);
        biomass = readgeoraster(bioname);
        biomis = bioinfo.MissingDataIndicator;
        if isempty(biomis)
            cl = bioinfo.NativeFormat;
            dum = mislist(ismember(cllist,cl));
            m1 = min(biomass,[],'all');
            m2 = max(biomass,[],'all');
            if m1==0, biomis = 0; end
            if m2==dum, biomis = [biomis,m2]; end %#ok<AGROW>
        end
        save(subfilename,'x','y','northlat','southlat','westlon','eastlon',...
            'biomass','bioinfo','biomis')
        clear biomass bioinfo ecomis
    else
        presbio(bb) = false; %#ok<SAGROW>
    end
    strcat("Done with importing ESA-CCI biomass data into block #",num2str(bb))
end

save(sparameterfile,'presbio','biomis','-append')


% Get the corresponding MODIS landcover
% *************************************
modisperblock = zeros(10,8,nbblocks);
modisfilesperblock = strings(nbblocks,10);
for bb = 1 : nbblocks
    if presbio(bb) == false, continue; end
    subfilename = strcat(subdata10,"BioSubData0008_",num2str(bb),".mat");
    load(subfilename,'northlat','southlat','westlon','eastlon','bioinfo')
    
    ll = logical(ismember(table2array(modiscoord(:,strcmp(modisvars,"lat_max"))),northlat) + ...
        ismember(table2array(modiscoord(:,strcmp(modisvars,"lat_min"))),southlat));
    kk = (table2array(modiscoord(:,strcmp(modisvars,"lon_min"))) <= westlon) + ...
        (table2array(modiscoord(:,strcmp(modisvars,"lon_max"))) >= westlon) == 2;
    nn = (table2array(modiscoord(:,strcmp(modisvars,"lon_max"))) >= eastlon) + ...
        (table2array(modiscoord(:,strcmp(modisvars,"lon_min"))) <= eastlon) == 2;
    
    oo = logical(kk+nn);
    mod = modiscoord(ll+oo==2,:); clear ll kk nn
    [nf,~] = size(mod);
    modisfiles = strings(nf,1);
    for kk = 1 : nf
        hp = num2str(table2array(mod(kk,strcmp(modisvars,"ih"))));
        if length(hp) == 1, hp = strcat('0',hp); end
        vp = num2str(table2array(mod(kk,strcmp(modisvars,"iv"))));
        if length(vp) == 1, vp = strcat('0',vp); end
        % check if the modif file actually exists
        modtile = strcat("h",hp,"v",vp);
        if sum(contains(modisfilelist0008,modtile)) == 1
            modisfiles(kk) = strcat("MLC10.",modtile,".tif");
        elseif sum(contains(modisoriginalfilelist,modtile)) == 1
            error(strcat('MODIS file with tile "',modtile,'" is missing from folder "',...
                modis0008folder,'"'))
        end
    end
    dd = modisfiles == ""; vv = logical(1-dd);
    modisfiles = modisfiles(vv); mod = mod(vv,:); nf = nf - sum(dd);
    modisfilesperblock(bb,1:nf) = modisfiles;
    % check "real" longitude limits as they might be different from original due to projection
    if nf > 0
        mfl = true(nf,1);
        lon = zeros(nf,2);
        lat = zeros(nf,2);
        for ff = 1 : nf
            mfname = strcat(modis0008folder,modisfiles(ff));
            modisinfo = georasterinfo(mfname);
            lon(ff,:) = modisinfo.RasterReference.LongitudeLimits;
            lat(ff,:) = modisinfo.RasterReference.LatitudeLimits;
            if lon(ff,2) <= westlon || lon(ff,1) >= eastlon
                mfl(ff) = false;
            end
            clear mfname modisinfo dum
        end
        modisfiles = modisfiles(mfl); nf = sum(mfl);
        if exist('mod','var') == 1
            mod = mod(mfl,:); lat = lat(mfl,:); lon = lon(mfl,:);
            modisperblock(1,1:4,bb) = [westlon,eastlon,southlat,northlat];
            modisperblock(2:nf+1,1:4,bb) = table2array(mod(:,3:6));
            modisperblock(2:nf+1,5:6,bb) = lon;
            modisperblock(2:nf+1,7:8,bb) = lat;
        end
        mfname = strcat(modis0008folder,modisfiles(1));
        modisinfo = georasterinfo(mfname);
        modisnodata = modisinfo.MissingDataIndicator;
        modisval = readgeoraster(mfname);
        if modisinfo.RasterSize(1) ~= bioinfo.RasterSize(1) || ...
                modisinfo.RasterSize(2) ~= bioinfo.RasterSize(2)
            modisval = snaparray(bioinfo,modisinfo,modisval);
        end
        if nf > 1
            for mm = 2 : nf
                mfname = strcat(modis0008folder,modisfiles(mm));
                modisinfo = georasterinfo(mfname);
                modisval2 = readgeoraster(mfname);
                if modisinfo.RasterSize(1) ~= bioinfo.RasterSize(1) || ...
                        modisinfo.RasterSize(2) ~= bioinfo.RasterSize(2)
                    modisval2 = snaparray(bioinfo,modisinfo,modisval2);
                end
                modisval(modisval==modisnodata) = modisval2(modisval==modisnodata);
                clear mfname modisinfo modisval2
            end
        end
        save(subfilename,'modisval','-append')
    else
        if northlat >= 0, lat = strcat(num2str(northlat),"N"); else, lat = strcat(num2str(abs(northlat)),"S"); end
        if westlon >= 0, lon = strcat(num2str(westlon),"E"); else, lon = strcat(num2str(abs(westlon)),"W"); end
        warning(strcat("there is no MODIS landcover values for block #",num2str(bb),...
            " [",lat,",",lon,"]"))
    end
    
    strcat("done with block #",num2str(bb)," - with ",num2str(nf)," MODIS files")
    clear subfilename northlat southlat westlon eastlon bioinfo ll kk nn oo mod nf modisfiles kk ...
        hp vp mfname modisinfo modisnodata modisval mm modisval2 dd
end
nlm = max(squeeze(sum(modisperblock(:,1,:)~=0)));
modisperblock08 = modisperblock(1:nlm,:,:);
modisfilesperblock08 = modisfilesperblock(:,1:nlm-1);
save(sparameterfile,'modisperblock08','modisfilesperblock08','-append')


% Get the corresponding ecoregion files (but loop through files instead)
% *************************************

% Get ecoregion coordinates
ecocoords = zeros(numel(ecofilelist),4);
for ee = 1 : numel(ecofilelist)
    ecofname = ecofilelist(ee);
    cn = char(extractBefore(ecofname,"_Ecoregions.tif"));
    if cn(1)=='N', northlat = str2double(cn(2:3)); else, northlat = -str2double(cn(2:3)); end
    if cn(4)=='E', westlon = str2double(cn(5:7)); else, westlon = -str2double(cn(5:7)); end
    southlat = northlat - 10;
    if numel(cn) > 7
        rem = cn(9:numel(cn));
        if rem(1) == 'E', eastlon = str2double(rem(2:4)); ...
        elseif cn(4)=='E', eastlon = str2double(rem(1:3)); else, eastlon = -str2double(rem(1:3)); end
        eastlon = eastlon + 10;
    else
        eastlon = westlon + 10;
    end
    % check that name corresponds to actual coordinates
    ecoinfo = georasterinfo(strcat(ecofolder,ecofname));
    nn = round(ecoinfo.RasterReference.LatitudeLimits);
    ww = round(ecoinfo.RasterReference.LongitudeLimits);
    if nn(1) == southlat && nn(2) == northlat && ww(1) == westlon && ww(2) == eastlon
        ecocoords(ee,:) = [northlat,southlat,westlon,eastlon];
    else
        warning(strcat('File "',ecofname,'" ',"has the wrong coordinates in it's name! ",...
            "Check that we have all necessary values in the end"))
        ecocoords(ee,:) = [nn(2),nn(1),ww];
    end
end
save(sparameterfile,"ecofilelist","ecocoords","-append")

% In this resolution, I seem to have used ObjectID as identifier
eco0008perblock = zeros(nbblocks,100);  nbmaxeco = 0;
for ee = 1 : numel(ecofilelist)
    ecofile = strcat(ecofolder,ecofilelist(ee));
    ecd = ecocoords(ee,:);
    sy = find(ncs == ecd(1));
    firstx = find(wcs == ecd(3));
    lastx = find(wcs == ecd(4)) - 1;
    xser = firstx : lastx;
    nxx = numel(xser);
    ecoinfo = georasterinfo(ecofile);
    ecomis = ecoinfo.MissingDataIndicator;
    ecoval = readgeoraster(ecofile);
    if isempty(ecomis)
        cl = ecoinfo.NativeFormat;
        ecomis = mislist(ismember(cllist,cl));
    end
    dum = unique(ecoval);
    if sum(~ismember(setdiff(dum,ecomis),ecoindex))> 0
        ecolist = setdiff(dum,ecomis);
        val = ecolist(~ismember(ecolist,ecoindex));
        vn = string(val);
        error(strcat("Check why I have value ",vn," in file ",ecofile))
    end
    dum = ecoinfo.RasterSize;
    if dum(2) == (blocksize0008 * nxx) + 1 % occasionally I have one too many pixels in longitude!
        lonlim = ecoinfo.RasterReference.LongitudeLimits;
        dd = lonlim - round(lonlim);
        side = find(dd == max(dd));
        if side == 1
            ecoval = ecoval(:,2:(blocksize0008*nxx+1));
        else
            ecoval = ecoval(:,1:blocksize0008*nxx);
        end
        dum = size(ecoval);
    end
    if dum(1) == blocksize0008 && dum(2)== blocksize0008 * nxx
        for xx = 1 : nxx
            sx = xser(xx);
            jj = (xx-1) * blocksize0008 + 1 : xx * blocksize0008;
            ecobio08 = ecoval(:,jj);
            bb = sx + (sy-1) * nblx;
            subfilename = strcat(subdata10,"BioSubData0008_",num2str(bb),".mat");
            load(subfilename,'x','y','biomass','biomis')
            if x == sx && y == sy
                % gap-fill if needed
                biomask = biomass~=biomis;
                a = find(ecobio08(biomask)==ecomis,1);
                if ~isempty(a)
                    ecobio08 = blockgapfill(ecobio08,biomask,blocksize0008,ecomis,50);
                    strcat("Gap-filled eco in block #",num2str(bb))
                end
                save(subfilename,'ecobio08','-append')
            else
                error("coordinate calculation is wrong")
            end
            bbecolist = setdiff(unique(ecobio08),ecomis);
            eco0008perblock(bb,1:numel(bbecolist)) = bbecolist;
            if numel(bbecolist) > nbmaxeco; nbmaxeco = numel(bbecolist); end
            clear sx jj ecobio08 bb subfilename x y bbecolist ej
        end
    else
        error("This should not have happened")
    end
    
    strcat("Done with ecofile ",ecofilelist(ee))
    clear ecofile ecoinfo ecoval ecd sy firstx lastx xser nxx dum xx
end
eco0008perblock = eco0008perblock(:,1:nbmaxeco);
save(sparameterfile,'eco0008perblock','-append')


% Calculate biomass statistics (loop through ecoregions)
% ****************************
ptsnames = strcat("ESABio",IGBPAbbBio(forwsaind));
ptsnames = cat(1,ptsnames,"ESABioAllForests");
esastats = zeros(numel(ecoindex),numel(ptsnames),numel(despct));
esaperecoid = zeros(numel(ecoindex),numel(ptsnames),numel(despct));
esapixperecoid = zeros(numel(ecoindex),numel(ptsnames));
% esastatsperbiome = zeros(nbiomes,numel(despct));
% esapixperbiome = zeros(nbiomes,1);

for bio = 1 : nbiomes
    tic
    dum = ecoperbiome(bio,:);
    thisbioecolist = dum(dum>0); clear dum
%     biopts = [];
    for ee = 1 : numel(thisbioecolist)
        thisecoid = thisbioecolist(ee);
        ecofile = strcat(subdata10,"Ecoregion#",num2str(thisecoid),"_0008.mat");
        thisecoindex = ecoindex(ecoid==thisecoid);
        ecoii = ecoindex == thisecoindex;
        [blocks,~] = find(ismember(eco0008perblock,thisecoindex));
        blocks = blocks(blocks>0);

        if numel(blocks) == 0, continue; end
        for ff = 1 : numel(ptsnames)
            eval(strcat(ptsnames(ff)," = [];"))
        end

        for bb = 1 : numel(blocks)
            bid = blocks(bb);
            if presbio(bid) == 0, continue, end
            subfilename = strcat(subdata10,"BioSubData0008_",num2str(bid),".mat");
            load(subfilename,'ecobio08','modisval','biomass','bioinfo')
            tilesize = bioinfo.RasterSize;
            kk = ecobio08 == thisecoindex;
            esa = nan(tilesize);
            esa(kk) = double(biomass(kk));
            esa(esa==0) = nan;
%             esa(esa==biomis) = nan; %looks like "missing values are set to zero
            esa = esa .* dm2C;
            for ff = 1 : numel(forwsaind) + 1
                if ff == numel(forwsaind) + 1, fi = forwsaind; else, fi = forwsaind(ff); end
                ii = ismember(modisval,fi);
                dum = esa(ii); vals = dum(isfinite(dum));
                eval(strcat(ptsnames(ff)," = cat(1,",ptsnames(ff),",vals);"))
                clear fi ii dum vals
            end
            clear ecobio08 modisval biomass bioinfo bid subfilename kk esa ii dum
        end
        for ff = 1 : numel(forwsaind) + 1
            eval(strcat("pts = ",ptsnames(ff),";"))
            stats = prctile(pts,despct);
            esastats(ecoii,ff,:) = stats;
            esaperecoid(thisecoid,ff,:) = stats;
            esapixperecoid(thisecoid,ff) = numel(pts);
            clear pts stats
        end
        save(ecofile,ptsnames{:})
        
%         biopts = cat(1,biopts,ESABioAllForests);
        clear(ptsnames{:})

        strcat("Done with statistics in ecoregion #",num2str(thisecoid),' "',econames(ecoii),'"')
        clear thisecoid thisecoindex blocks ff bb
    end
    
%     esastatsperbiome(bio,:) = prctile(biopts,despct);     % too much memory!
%     esapixperbiome(bio) = numel(biopts);
    clear thisbioecolist biopts ee
    strcat("Done with statistics in biome -",biomenames(bio))
    eltime = toc;
    if eltime > 180
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
end

% save(ReforOppfname,'esastats','esaperecoid','esapixperecoid','esastatsperbiome','esapixperbiome','-append')
save(ReforOppfname,'esastats','esaperecoid','esapixperecoid','-append')