% Global Albedo - Kernel Statistics
%
% Uses results from the GlobalAlbedo Glasgow COP runs, namely outputs from different kernels data
%   and compiles some statistics. It also computes statistics on the original input kernels data.
%
% Created 10/19/2021 by Natalia Hasler, Clark University

% Last modified: NH 10/19/2021


clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A. Define and load variables
% *****************************


% **********************
% 1. User-Defines variables
% -------------------------

computer = "GEO-007264";                % where the routine is run (for directory paths):
                                        %   GEO-007264 (my laptop)
                                        %   GEO-005199 (old supercomputer)
                                        %   GEO-006381 (new supercomputer)
desiredstatistics = ...                 % list of statistics to be performed
    ["min","max","mean","median","range","standard deviation"];

% **********************


% 2. Define inpouts
% -----------------
switch computer
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
        mapdir = strcat('G:\GlobalAlbedo\');
        resultdir = rootdir;
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
        mapdir = strcat('G:\TNC\GlobalAlbedo\');
        resultdir = 'E:\NCS-GlobalAlbedo\';
end

glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'IGBPBiomes','biomes_igbp','latlonscale','nblocks','nlat','nlon','nlcc',...
    'nmonths','nkernels','npw','presland','pathwayslandcover','pwnames','pwid','kernels',...
    'regoutputfiles','resultslocalfolder','landmask','blocksize')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% B. Compute TOA RF and CO2e of land conversion scenarios
% ********************************************************

% 1. Loop through sub-regions
% ---------------------------
tilediv = 4;
tlsize = blocksize / tilediv;
nbtiles = tilediv * tilediv;
        
for rr = 1 : nblocks
    if presland(rr) == true
        subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
        load(subdatafname);
        
        
        % 2. Sub-divide region in smaller tiles for use in parfor loop
        % ------------------------------------------------------------
        tilemvars = NaN(tlsize,tlsize,nmonths,nbtiles);
        tilealvar = NaN(tlsize,tlsize,nmonths,nbiomes,nbtiles);
        tileco2e = NaN(tlsize,tlsize,nkernels,nlcc,nbtiles);
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

        
        % 3. Loop through tiles, define arrays and find land pixels
        % ---------------------------------------------------------
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
            sm_pathwaymiss = false(tlsize,tlsize,nmonths,nlcc);
            
            [a,b] = find(sm_landmask);
        
        
            % 4. Loop through land pixels/pathways and calculate monthly radiative forcings
            % -----------------------------------------------------------------------------
            for pt = 1 : length(a)
                i = a(pt);
                j = b(pt);
                
                % a. Loop through individual pathways
                for ll = 1 : nlcc
                    lc1 = biomes_igbp == pathwayslandcover(ll,1); %#ok<PFBNS>
                    lc2 = biomes_igbp == pathwayslandcover(ll,2);
                    ta = zeros(nmonths,nkernels);
                    
                    % b. Loop through months                               
                    for m = 1 : nmonths
                    
                        % c. Calculate surface albedo differences
                        daNSBS = sm_nosnow_blacksky(i,j,m,lc1) - sm_nosnow_blacksky(i,j,m,lc2);
                        daNSWS = sm_nosnow_whitesky(i,j,m,lc1) - sm_nosnow_whitesky(i,j,m,lc2);
                        daSCBS = sm_snowcov_blacksky(i,j,m,lc1) - sm_snowcov_blacksky(i,j,m,lc2);
                        daSCWS = sm_snowcov_whitesky(i,j,m,lc1) - sm_snowcov_whitesky(i,j,m,lc2);

                        % d. Calculate monthly Surface Weighted Albedo Differences
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
                                        
                        % e. Compute monthly TOA RF in [W m-2]
                        for kk = 1 : nkernels
                            ta(m,kk) =  Wda .* sm_kernels(i,j,m,kk) .*kernelscale(kk); %#ok<PFBNS>
                        end
                    
                    end     % end of months loop
                
                
                    % 5. Save missing values for quality control/improvements
                    % -------------------------------------------------------
                    %    (Kernels don't have missing values, so I can ignore that dimension)
                    sm_pathwaymiss(i,j,:,ll) = isnan(squeeze(ta(:,1)));
                
                
                    % 6. Calculate CO2 equivalent
                    % ---------------------------
                    % a. Annual mean RF
                    arf = squeeze(mean(ta,1,'omitnan'));
                    
                    % b. Global annual mean RF
                    glob = arf .* sm_pixelarea(i,j) ./ Aglobe;
                
                    % c. CO2e (Pg C)
                    %    (also, flip sign for ease of read)
                    PgC_TOA = Ca_PgC .* exp(-1 .* glob ./ 5.35);
                                    
                    % d. in tonnes CO2 ha-2
                    CO2eq_TOA = (PgC_TOA - Ca_PgC) .* 1e13 ./ sm_pixelarea(i,j) .* 44/12;
                
                    % e. Save values in tile arrays
                    sm_pathwayco2(i,j,:,ll) = CO2eq_TOA;
                
                end     % end loop on all individual pathways
                                
            end     % end loop on all forested pixels
            
            
            % 8. Save tile results back into regional arrays
            % ----------------------------------------------
            tileco2e(:,:,:,:,tl) = sm_pathwayco2;
            tilemiss(:,:,:,:,tl) = sm_pathwaymiss;
            
            strcat("done with RF/CO2e calculation for subregion ",num2str(rr)," - tile #",...
                num2str(tl)," (",num2str(length(a))," points)")

        end     % end parfor on tiles
        
        
        % 9. Recombine tiles into one regional arrays
        % -------------------------------------------        
        for pp = 1 : npw
            if sum(pwid==pp) == nfrst
                eval(strcat("regco2",lower(pwnames(pp))," = regfCO2e;"));
                eval(strcat("regmiss",lower(pwnames(pp))," = regfmiss;"));
            else
                eval(strcat("regco2",lower(pwnames(pp))," = regnfCO2e;"));
                eval(strcat("regmiss",lower(pwnames(pp))," = regnfmiss;"));
            end
        end
        
        if ismember(rr,canblocks)
            canrf = regcanmonths; canda = regcanmonths; cangra = regcan; cancro = regcan;
        end

        for tl = 1 : nbtiles
            y = ceil(tl/tilediv);
            x = tl - (y-1)*tilediv;
            tlj = (x-1)*tlsize + 1 : x * tlsize;
            tli = (y-1)*tlsize + 1 : y * tlsize;
            
            for pp = 1 : npw
                ii = pwid == pp;
                eval(strcat("regco2",lower(pwnames(pp)),"(tli,tlj,:,:) = tileco2e(:,:,:,ii,tl);"));
                eval(strcat("regmiss",lower(pwnames(pp)),"(tli,tlj,:,:) = tilemiss(:,:,:,ii,tl);"));
            end
        end            
        
        
        % 10. Save values in global arrays
        % --------------------------------
        for pp = 1 : npw
            eval(strcat(pwnames(pp),"(regi,regj,:,:) = regco2",lower(pwnames(pp)),";"));
            eval(strcat("Miss",pwnames(pp),"(regi,regj,:,:) = regmiss",lower(pwnames(pp)),";"));
        end

        regvars = cellstr(cat(1,strcat("regco2",lower(pwnames)),...
            strcat("regmiss",lower(pwnames))));
        save(subdatafname,regvars{:},'-append')
        
        strcat("done with RF/CO2e calculation for subregion ",num2str(rr),...
            " (",num2str(sum(blocklandmask,'all'))," points)")
        
        clear(regvars{:})
        clear nosnowblacksky nosnowwhitesky snowcovblacksky snowcovwhitesky ...
            x y regi regj landcoverprop blocklandmask snowcover directsolar ...
            diffusesolar diffusefraction cam3 cam5 echam6 hadgem2 ...
            canrf canda cangra cancro tlnsbs tlnsws tlscbs tlscws tlsnowcover tldiffrac ...
            tlkernels tllandmask tlpixarea tileco2e tilemiss ...
            tlcanrf tlcanad tlcangra tlcancro

        
    end     % end if region has land data
end
        
datafname = strcat(regoutputfiles,"AlbedoFinalData6ker.mat");
globalvars = cellstr(cat(1,pwnames,strcat("Miss",pwnames)));
save(datafname,globalvars{:},'-v7.3')

