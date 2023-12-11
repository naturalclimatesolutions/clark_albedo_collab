% Albedo-change induced TOA Radiative Forcing for single land-cover pathway
%   (if all land was one land-cover changing to a single different land-cover)


% Define variables
% ****************
pathways = ...                          % name of studied pathways
    ["deforestation to croplands";"deforestation to grasslands";...
    "deforestation to urban developement";"reforestation from open shrublands";...
    "reforestation from closed shrublands";"crop mosaic to forest"];
pwnames = ["Forest2Crop";...            % pathway names to use as folder names
    "Forest2Grass";"Forest2Urban";...
    "OpenShrub2Forest";"ClosedShrub2Forest";"CropVegMosaic2Forest"];
npw = length(pathways);
nfrst = length(forestsavannabiomes);


plc = strings(nfrst,3,npw);           % String array of 3-letter acronyms of old/new land cover in
                                        %   the desired pathway (one line) with 3 columns as
                                        %   1.  3-letter acronym of initial land cover
                                        %   2.  3-letter acronym of final land cover
                                        %   3.  pathway number (matching line number in "pathways")
plc(:,1,1:3) = repmat(forestsavannabiomes,[1,1,3]);
plc(:,2,1) = repmat("CRO",[nfrst,1]);
plc(:,2,2) = repmat("GRA",[nfrst,1]);
plc(:,2,3) = repmat("URB",[nfrst,1]);
plc(:,2,4:6) = repmat(forestsavannabiomes',[1,1,3]);
plc(:,1,4) = repmat("OSH",[nfrst,1]);
plc(:,1,5) = repmat("CSH",[nfrst,1]);
plc(:,1,6) = repmat("MOS",[nfrst,1]);
plc(:,3,:) = permute(repmat(1:6,[nfrst,1]),[1,3,2]);
pathwayslandcover = reshape(permute(plc,[1,3,2]),[nfrst*(npw),3]);
nlcc = length(pathwayslandcover);       % number of individual land cover conversions



% Calculate Radiative Forcing
% ***************************
tilediv = 4;
tlsize = blocksize05 / tilediv;
nbtiles = tilediv * tilediv;
        
for rr = 1 : nblocks
    if presland(rr) == false, continue; end
    tic
    subdatafname = strcat(regoutputfiles,"Albedo05_",num2str(rr),".mat");
    load(subdatafname);


    % Sub-divide region in smaller tiles for use in parfor loop
    tilemvars = NaN(tlsize,tlsize,nmonths,nbtiles);
    tilealvar = NaN(tlsize,tlsize,nmonths,numel(IGBPBiomes),nbtiles);
    tileco2e = NaN(tlsize,tlsize,nkernels,nlcc,nbtiles);
    tileco2estat = nan(tlsize,tlsize,nstats,nlcc,nbtiles);
    tilemiss = false(tlsize,tlsize,nmonths,nlcc,nbtiles);

    tlnsbs = tilealvar; tlnsws = tilealvar; tlscbs = tilealvar; tlscws = tilealvar;
    tlsnowcover = tilemvars; tldiffrac = tilemvars;

    tlkernels = NaN(tlsize,tlsize,nmonths,nkernels,nbtiles);
    tllandmask = false(tlsize,tlsize,nbtiles);
    tlpixarea = NaN(tlsize,tlsize,nbtiles);

    clear tilemvars tilealvar

    for tl = 1 : nbtiles
        y = ceil(tl/tilediv);
        x = tl - (y-1)*tilediv;
        tlj = (x-1)*tlsize + 1 : x * tlsize;
        tli = (y-1)*tlsize + 1 : y * tlsize;

        tlnsbs(:,:,:,:,tl) = nosnowblacksky(tli,tlj,:,:);
        tlnsws(:,:,:,:,tl) = nosnowwhitesky(tli,tlj,:,:);
        tlscbs(:,:,:,:,tl) = snowcovblacksky(tli,tlj,:,:);
        tlscws(:,:,:,:,tl) = snowcovwhitesky(tli,tlj,:,:);

        tlsnowcover(:,:,:,tl) = snowcover(tli,tlj,:);
        tldiffrac(:,:,:,tl) = diffusefraction(tli,tlj,:);

        tllandmask(:,:,tl) = blocklandmask(tli,tlj);
        tlpixarea(:,:,tl) = pixelarea(tli,tlj);

        for kk = 1 : nkernels
            eval(strcat("kern = ",lower(kernels(kk)),";"))
            tlkernels(:,:,:,kk,tl) = kern(tli,tlj,:);
        end

        clear x y tli tlj kk
    end


    % Loop through tiles, define arrays and find land pixels
    parfor tl = 1 : nbtiles

        sm_nosnow_blacksky = tlnsbs(:,:,:,:,tl);
        sm_nosnow_whitesky = tlnsws(:,:,:,:,tl);
        sm_snowcov_blacksky = tlscbs(:,:,:,:,tl);
        sm_snowcov_whitesky = tlscws(:,:,:,:,tl);

        sm_snowcover = tlsnowcover(:,:,:,tl);
        sm_diffuse_fraction = tldiffrac(:,:,:,tl);
        sm_kernels = tlkernels(:,:,:,:,tl);

        sm_landmask = tllandmask(:,:,tl);
        sm_pixelarea = tlpixarea(:,:,tl);

        sm_pathwayco2 = NaN(tlsize,tlsize,nkernels,nlcc);
        sm_pathwayco2stat = nan(tlsize,tlsize,nstats,nlcc);
        sm_pathwaymiss = false(tlsize,tlsize,nmonths,nlcc);

        [a,b] = find(sm_landmask);


        % Loop through land pixels/pathways and calculate monthly radiative forcings
        for pt = 1 : length(a)
            i = a(pt);
            j = b(pt);

            % Loop through individual pathways
            for ll = 1 : nlcc
                lc1 = IGBPAbbBio == pathwayslandcover(ll,1); %#ok<PFBNS>
                lc2 = IGBPAbbBio == pathwayslandcover(ll,2);
                ta = zeros(nmonths,nkernels);

                % Loop through months
                for m = 1 : nmonths

                    % Calculate surface albedo differences
                    daNSBS = sm_nosnow_blacksky(i,j,m,lc1) - sm_nosnow_blacksky(i,j,m,lc2);
                    daNSWS = sm_nosnow_whitesky(i,j,m,lc1) - sm_nosnow_whitesky(i,j,m,lc2);
                    daSCBS = sm_snowcov_blacksky(i,j,m,lc1) - sm_snowcov_blacksky(i,j,m,lc2);
                    daSCWS = sm_snowcov_whitesky(i,j,m,lc1) - sm_snowcov_whitesky(i,j,m,lc2);

                    % Calculate monthly Surface Weighted Albedo Differences
                    %    (to avoid NaN values due to undefined albedos in irrelevant situations,
                    %     I add a cumbersome if statement ... now that we filled these, it
                    %     should not be necessary, but I left it as is)
                    snow = sm_snowcover(i,j,m);
                    if snow == 0
                        Wda = squeeze(daNSBS .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                            daNSWS .* sm_diffuse_fraction(i,j,m));
                    elseif snow == 1
                        Wda = squeeze(daSCBS .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                            daSCWS .* sm_diffuse_fraction(i,j,m));
                    else
                        Wda = squeeze(...
                            daNSBS .* (1 - snow) .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                            daNSWS .* (1 - snow) .* sm_diffuse_fraction(i,j,m) + ...
                            daSCBS .* snow .* (1 - sm_diffuse_fraction(i,j,m)) + ...
                            daSCWS .* snow .* sm_diffuse_fraction(i,j,m));
                    end

                    % Compute monthly TOA RF in [W m-2]
                    for kk = 1 : nkernels
                        ta(m,kk) =  Wda .* sm_kernels(i,j,m,kk) .*kernelscale(kk);
                    end

                end     % end of months loop


                sm_pathwaymiss(i,j,:,ll) = isnan(squeeze(ta(:,1)));


                % Calculate CO2 equivalent
                arf = squeeze(mean(ta,1,'omitnan'));        % Annual mean RF
                glob = arf .* sm_pixelarea(i,j) ./ Aglobe;  % Global annual mean RF
                PgC_TOA = Ca_PgC .* exp(-1 .* glob ./ 5.35);% CO2e (Pg C)
                CO2eq_TOA = (PgC_TOA - Ca_PgC) .* 1e13 ...  % in tonnes CO2 ha-2
                    ./ sm_pixelarea(i,j) .* 44/12;
                CO2eqTOAGWP = CO2eq_TOA .* GWP100;          % accounting for land/ocean response
                sm_pathwayco2(i,j,:,ll) = CO2eqTOAGWP;


                % Calculate kernel's statistics
                % min/max values should represent the min and max offset, not the minimum value,
                % so it gets tricky when values cross the line, but I'll keep it somewhat simple
                % with generally max effect meaning minimum value.
                vsig = sign(CO2eqTOAGWP);
                if sum(vsig) == 6
                    minval = min(CO2eqTOAGWP); maxval = max(CO2eqTOAGWP);
                else
                    maxval = min(CO2eqTOAGWP);
                    mi = min(abs(CO2eqTOAGWP));
                    ii = abs(CO2eqTOAGWP) == mi;
                    tmpminval = CO2eqTOAGWP(ii);
                    if numel(tmpminval) == 1
                        minval = tmpminval;
                    elseif sum(diff(tmpminval)) == 0
                        minval = tmpminval(1);
                    else
                        minval = mi;
                    end
                end
                for ss = 1 : nstats
                    stat = singlepathstatistics(ss);
                    switch stat
                        case "min"
                            sm_pathwayco2stat(i,j,ss,ll) = minval;
                        case "max"
                            sm_pathwayco2stat(i,j,ss,ll) = maxval;
                        case "median"
                            sm_pathwayco2stat(i,j,ss,ll) = median(CO2eqTOAGWP,'omitnan');
                        case "mean"
                            sm_pathwayco2stat(i,j,ss,ll) = mean(CO2eqTOAGWP,'omitnan');
                        case "standard deviation"
                            sm_pathwayco2stat(i,j,ss,ll) = std(CO2eqTOAGWP,'omitnan');
                        case "range"
                            sm_pathwayco2stat(i,j,ss,ll) = range(CO2eqTOAGWP);
                    end
                end

            end     % end loop on all individual pathways
        end     % end loop on all forested pixels


        % Save tile results back into regional arrays
        tileco2e(:,:,:,:,tl) = sm_pathwayco2;
        tilemiss(:,:,:,:,tl) = sm_pathwaymiss;
        tileco2estat(:,:,:,:,tl) = sm_pathwayco2stat;

%         strcat("done with block #",num2str(rr)," - tile #",...
%             num2str(tl)," (",num2str(length(a))," points)")

    end     % end parfor on tiles


    % Recombine tiles into one regional arrays
    regco2e = nan(blocksize05,blocksize05,nkernels,nlcc);
    regco2estat = nan(blocksize05,blocksize05,nstats,nlcc);
    regmiss = false(blocksize05,blocksize05,nmonths,nlcc);

    for tl = 1 : nbtiles
        y = ceil(tl/tilediv);
        x = tl - (y-1)*tilediv;
        tlj = (x-1)*tlsize + 1 : x * tlsize;
        tli = (y-1)*tlsize + 1 : y * tlsize;

        regco2e(tli,tlj,:,:) = tileco2e(:,:,:,:,tl);
        regco2estat(tli,tlj,:,:) = tileco2estat(:,:,:,:,tl);
        regmiss(tli,tlj,:,:) = tilemiss(:,:,:,:,tl);
    end

    save(subdatafname,"regco2e","regmiss","regco2estat",'-append')

    clear nosnowblacksky nosnowwhitesky snowcovblacksky snowcovwhitesky ...
        x y regi regj landcoverprop snowcover directsolar ...
        diffusesolar diffusefraction cam3 cam5 echam6 hadgem2 ...
        canrf canda cangra cancro tlnsbs tlnsws tlscbs tlscws tlsnowcover tldiffrac ...
        tlkernels tllandmask tlpixarea tileco2e tilemiss ...
        tlcanrf tlcanad tlcangra tlcancro regco2e regmiss regco2estat

    strcat("done with RF/CO2e calculation for subregion ",num2str(rr),...
        " (",num2str(sum(blocklandmask,'all'))," points)")
    elmtime = toc;
    if elmtime > 180
        h = floor(elmtime/3600);
        m = floor((elmtime - h*3600)/60);
        s = elmtime - (h*3600 + m*60);
        duration([h,m,s])
    end
    clear elmtime h m s blocklandmask
end