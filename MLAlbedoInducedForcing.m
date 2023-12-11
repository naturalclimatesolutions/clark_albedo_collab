% Albedo-Induced Climate Forcing (from "most likely open lands" to "most likely forests"

% Prepare variables
% *****************
kervars = lower(kernels);
othvars = ["snowcover","diffusefraction"];
inputvars = [kervars,othvars,albsn]';
outputvars = strcat("tiled",inputvars);
bigtosmall = round(biggestblock./bigblocks); % round to avoid warnings
bigblocksize005 = biggestblock ./ latlon005;
downscaletilesize = latlon05 ./ latlon005;
npfx = 5; npfy = 5; npftiles = npfx * npfy;
pftilesize = blocksize005 ./ npfx;
bigtosmallblockscorrespondance = zeros(nbblocks,4);
kerstat = ["median","min","max"];
% NOTE : if kerstat changes than change the entire code accordingly!!!

nigbp = numel(LC005ind);
for rr = 1 : nblocks
    y = ceil(rr/nblx);
    x = rr - (y-1)*nblx;
    x2 = (x-1)*bigtosmall + 1 : x*bigtosmall;
    y2 = (y-1)*bigtosmall + 1 : y*bigtosmall;
    for xx = 1 : bigtosmall
        sx = x2(xx);
        for yy = 1 : bigtosmall
            sy = y2(yy);
            bb = sx + (sy-1) * nnblx;
            bigtosmallblockscorrespondance(bb,:) = [rr,bb,yy,xx];
        end
    end
end

% Downscale the base data to the finer grid (and tile it)
% *****************************************
for rr = 1 : nblocks
    tic
    if presland(rr) == false, continue; end
    subdatafname = strcat(regoutputfiles,"Albedo05_",num2str(rr),".mat");
    load(subdatafname,inputvars{:});

    bli = bigtosmallblockscorrespondance(:,1) == rr;
    blocks = bigtosmallblockscorrespondance(bli,2);
    smx = bigtosmallblockscorrespondance(bli,4);
    smy = bigtosmallblockscorrespondance(bli,3);


    for vv = 1 : numel(inputvars)
        eval(strcat("var = ",inputvars(vv),";"))
        [~,~,~,nlc] = size(var);

        if nlc == nigbp
            var = var(:,:,:,[2:nigbp,1]);       % albedo has water as the first index
            newvar = uint16(var * 1000);
            var005 = zeros(bigblocksize005,bigblocksize005,nmonths,nlc,'uint16');
            tiledvar = zeros(pftilesize,pftilesize,nmonths,nlc,npftiles,'uint16');

        elseif nlc == 1
            newvar = var;
            var005 = zeros(bigblocksize005,bigblocksize005,nmonths);
            tiledvar = zeros(pftilesize,pftilesize,nmonths,nlc,npftiles);
        else
            error(strcat("I have an unexpected dimension in variable ",inputvars(vv)))
        end
        for j = 1 : blocksize05
            jj = (j-1) * downscaletilesize + 1 : j*downscaletilesize;
            for i = 1 : blocksize05
                ii = (i-1) * downscaletilesize + 1 : i*downscaletilesize;
                var005(ii,jj,:,:) = repmat(newvar(i,j,:,:),[downscaletilesize,downscaletilesize,1]);
            end
        end

        for z = 1 : numel(blocks)
            bid = blocks(z);
            if preswlk(bid) == false, continue; end
            roinfname = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
            load(roinfname,"landmask")
            alfname = strcat(regoutputfiles,"Albedo005_",num2str(bid),".mat");
            ii = (smy(z)-1) * blocksize005 + 1 : smy(z)*blocksize005;
            jj = (smx(z)-1) * blocksize005 + 1 : smx(z)*blocksize005;
            smallvar = var005(ii,jj,:,:);
            if strcmp(inputvars(vv),"snowcover")
                smallvar = blockgapfill(smallvar,landmask,blocksize005,nan,50);  % gap-fill snow
            end
            for tl = 1 : npftiles
                y = ceil(tl/npfx);
                x = tl - (y-1)*npfx;
                j = (x-1)*pftilesize + 1 : x*pftilesize;
                i = (y-1)*pftilesize + 1 : y*pftilesize;
                tiledvar(:,:,:,:,tl) = smallvar(i,j,:,:);
            end
            ntiledvar = squeeze(tiledvar);
            eval(strcat(inputvars(vv)," = smallvar;"))
            eval(strcat(outputvars(vv)," = ntiledvar;"))
            if exist(alfname,'file') == 0
                eval(strcat("save(alfname,'",inputvars(vv),"','",outputvars(vv),"')"))
            else
                eval(strcat("save(alfname,'",inputvars(vv),"','",outputvars(vv),"','-append')"))
            end
        end
    end

    clear(inputvars{:},outputvars{:})
    clear var nlc newvar tiledvar i j ii jj z bid alfname smallvar x y ntiledvar blocks bli smy smx

    strcat("done with downscaling albedo variables in big block #", num2str(rr))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end

end


% Calculate radiative forcing (parfor through tiles)
% ***************************
for bb = 1 : nbblocks
    tic
    if preswlk(bb) == false, continue; end
    mlfname = strcat(regoutputfiles,"MostLikely005_",num2str(bb),".mat");
    load(mlfname,"MLOL","MLFC")
    rofname = strcat(regoutputfiles,"ROinputs_",num2str(bb),".mat");
    load(rofname,'pixarea','landmask')

    tiledMLOL = ones(pftilesize,pftilesize,npftiles,'uint8') .* u8noval;
    tiledMLFC = ones(pftilesize,pftilesize,npftiles,'uint8') .* u8noval;
    tiledpixarea = zeros(pftilesize,pftilesize,npftiles);
    tiledlandmask = false(pftilesize,pftilesize,npftiles);
    datatiles = true(npftiles,1);

    for tl = 1 : npftiles
        y = ceil(tl/npfx);
        x = tl - (y-1)*npfx;
        j = (x-1)*pftilesize + 1 : x*pftilesize;
        i = (y-1)*pftilesize + 1 : y*pftilesize;
        tiledMLOL(:,:,tl) = MLOL(i,j);
        tiledMLFC(:,:,tl) = MLFC(i,j);
        tiledpixarea(:,:,tl) = pixarea(i,j);
        mask = landmask(i,j);
        if sum(mask,'all') == 0, datatiles(tl) = false; end
        tiledlandmask(:,:,tl) = mask;
    end
    
    clear MLOL MLFC pixarea landmask

    alfname = strcat(regoutputfiles,"Albedo005_",num2str(bb),".mat");
    load(alfname,outputvars{:},"snowcover")
    for tl = 1 : npftiles % Because I gap-filled snowcover after the fact
        y = ceil(tl/npfx);
        x = tl - (y-1)*npfx;
        j = (x-1)*pftilesize + 1 : x*pftilesize;
        i = (y-1)*pftilesize + 1 : y*pftilesize;
        tiledsnowcover(:,:,:,tl) = snowcover(i,j,:); %#ok<SAGROW>
    end

    tiledkernels = zeros(pftilesize,pftilesize,nmonths,nkernels,npftiles);
    for kk = 1 : nkernels
        eval(strcat("tiledkernels(:,:,:,kk,:) = permute(tiled",kervars(kk),...
            ",[1:3,5,4]) .* kernelscale(kk);"))
    end
    clear tiledcam3 tiledcam5 tiledhadgem3 tiledhadgem2 tiledecham6 tiledcack1

    RFmed005 = nan(blocksize005,blocksize005);
    RFmin005 = nan(blocksize005,blocksize005);
    RFmax005 = nan(blocksize005,blocksize005);

    tiledRFmed = nan(pftilesize,pftilesize,npftiles);
    tiledRFmin = nan(pftilesize,pftilesize,npftiles);
    tiledRFmax = nan(pftilesize,pftilesize,npftiles);

    parfor tl = 1 : npftiles
        if datatiles(tl) == false, continue; end
        tlRFmed = tiledRFmed(:,:,tl);
        tlRFmin = tiledRFmin(:,:,tl);
        tlRFmax = tiledRFmax(:,:,tl);
        tlmlol = tiledMLOL(:,:,tl);
        tlmlfc = tiledMLFC(:,:,tl);
        tlpix = tiledpixarea(:,:,tl);
        tlmask = tiledlandmask(:,:,tl);
        tlsnow = squeeze(tiledsnowcover(:,:,:,tl));
        tldiff = squeeze(tileddiffusefraction(:,:,:,tl));
        tlnosnowblsk = squeeze(tilednosnowblacksky(:,:,:,:,tl));
        tlnosnowwhsk = squeeze(tilednosnowwhitesky(:,:,:,:,tl));
        tlsnowcovblsk = squeeze(tiledsnowcovblacksky(:,:,:,:,tl));
        tlsnowcovwhsk = squeeze(tiledsnowcovwhitesky(:,:,:,:,tl));
        tlkern = squeeze(tiledkernels(:,:,:,:,tl));

        [a,b] = find(tlmask);
        
        for px = 1 : numel(a)
            i = a(px);
            j = b(px);
            lc1 = tlmlol(i,j);
            lc2 = tlmlfc(i,j);
            if ~ismember(lc1,openlandind) || ~ismember(lc2,forwsaind), continue; end
            area = tlpix(i,j);
            nsbs = double(squeeze(tlnosnowblsk(i,j,:,:))) ./ 1000;
            nsws = double(squeeze(tlnosnowwhsk(i,j,:,:))) ./ 1000;
            scbs = double(squeeze(tlsnowcovblsk(i,j,:,:))) ./ 1000;
            scws = double(squeeze(tlsnowcovwhsk(i,j,:,:))) ./ 1000;
            snow = squeeze(tlsnow(i,j,:));
            diffuse = squeeze(tldiff(i,j,:));
            kerval = squeeze(tlkern(i,j,:,:));
            
            pval = ptradforcing(lc1,lc2,nsbs,nsws,scbs,scws,snow,diffuse,kerval,area,...
                kerstat,GWP100,Ca_PgC);

            tlRFmed(i,j) = pval(1);     % According to "kerstat"
            tlRFmin(i,j) = pval(2);
            tlRFmax(i,j) = pval(3);
        end
        
        tiledRFmed(:,:,tl) = tlRFmed;
        tiledRFmin(:,:,tl) = tlRFmin;
        tiledRFmax(:,:,tl) = tlRFmax;
        
    end
    
    % Save Radiative forcing
    % **********************
    for tl = 1 : npftiles
        y = ceil(tl/npfx);
        x = tl - (y-1)*npfx;
        j = (x-1)*pftilesize + 1 : x*pftilesize;
        i = (y-1)*pftilesize + 1 : y*pftilesize;
        RFmed005(i,j) = tiledRFmed(:,:,tl);
        RFmin005(i,j) = tiledRFmin(:,:,tl);
        RFmax005(i,j) = tiledRFmax(:,:,tl);
    end

    GWPmed005 = RFmed005 .* GWP100;
    GWPmin005 = RFmin005 .* GWP100;
    GWPmax005 = RFmax005 .* GWP100;
    
    
    radforfname = strcat(regoutputfiles,"RadForcing005_",num2str(bb),".mat");
    save(radforfname,"RFmed005","RFmin005","RFmax005","GWPmax005","GWPmin005","GWPmed005");
    
    
    clear RFmed005 RFmin005 RFmax005 rofname radforfname z bid ii jj  tl y x i j ...
        tiledkernels tiledRFmax tiledRFmin tiledRFmed
    clear(inputvars{:},outputvars{:})

    strcat("Done with Albedo-induced radiative forcing calculation in block #",num2str(bb))
    eltime = toc;
    if eltime > 120
        h = floor(eltime/3600);
        m = floor((eltime  - h*3600)/60);
        s = eltime - (h*3600 + m*60);
        duration([h,m,s])
    end
    
end
RFmaps = ["RFmed005","RFmin005","RFmax005"];
GWPmaps = ["GWPmed005","GWPmin005","GWPmax005"];
save(sparameterfile,"RFmaps","GWPmaps",'-append')






