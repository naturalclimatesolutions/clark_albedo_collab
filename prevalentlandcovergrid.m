function prevalentlc = prevalentlandcovergrid(resolution,landcover,ecoregion,pixelarea,...
    landmask,landcoverindex,mostlikelypereco,prevalentlandcoverthreshold,...
    combinesavannaflag,latlonscale,econoval)

% FUNCTION prevalentlandcovergrid: find prevalent landcover type in a grid cell
%   For each grid cell, this function determines the prevalent landcover-type of a predetermined
%   subset of landcover types (landcoverindex) in all land pixels within the grid cell that belong
%   to the same ecoregion.
%
% Input variables are as follows:
%   1.  resolution                  Resolution of the "grid" map or grid cell resolution
%   2.  landcover                   modis land cover map(s) - 
%                                       if several, saved in structure array in priority order
%   3.  ecoregion                   ecoregion map
%   4.  pixelarea                   map of orignal pixel areas
%   5.  landmask                    land pixels (non-water)
%   6.  landcoverindex              index list of landcover of interest
%   7.  mostlikelypereco            most likely land cover type per ecoregion
%   8.  prevalentlandcoverthreshold threshold below which presence of specific land cover is
%                                   ignored. Default is 1%
%   9.  combinesavannafalg          If true, combine both savanna land cover into one.
%  10.  latlonscale                 resolution of the original maps. Default is 0.005 decimal degrees
%  11.  econoval                    no-data value for ecoregion map. Default is 65535
%
% Note: I now constrain the prevalent land cover to areas that are at least > 1% or 20 pixels
%       (whichever is smaller) of the target area


u8noval = 2^8-1;
if nargin < 11
    econoval = 2^16-1;
end
if nargin < 10
    latlonscale = 0.005;
end
if nargin < 9
    combinesavannaflag = false;
end
if nargin < 8
    prevalentlandcoverthreshold = 1;
end

[blocksize,~] = size(ecoregion);

if isa(landcover,'struct')
    modisvarnames = string(fields(landcover));
    nlcm = numel(modisvarnames);
    dum =  ones(blocksize,blocksize,nlcm,'uint8') .* u8noval;
    for yy = 1 : nlcm
        eval(strcat(modisvarnames(yy)," = landcover.",modisvarnames(yy),";"))
        eval(strcat("dum(:,:,yy) = ",modisvarnames(yy),";"))
    end
    landcover = dum;
else
    modisvarnames = "landcover";
end


tilesize = resolution / latlonscale;
nblx = blocksize / tilesize;
nbly = blocksize / tilesize;
ntiles = nblx * nbly;

prevalentlc = ones(blocksize,blocksize,'uint8') .* u8noval;
prevalentlc(landmask) = 0;

if combinesavannaflag == true
    landcover(landcover==9) = 8;
end


if ntiles > 100
    nclx = 10;
    clustersize = blocksize / nclx;
    if ceil(clustersize/tilesize) == clustersize/tilesize
        ntlx = clustersize/tilesize;
        ntlpcl = ntlx * ntlx;
        nclusters = nclx * nclx;
    else
        error("code this")
    end
    
    cllandcover = ones(clustersize,clustersize,numel(modisvarnames),nclusters,'uint8') .* u8noval;
    clprevalentlc = ones(clustersize,clustersize,nclusters) .* u8noval;
    clecoregion = zeros(clustersize,clustersize,nclusters);
    clpixarea = zeros(clustersize,clustersize,nclusters);
    cllandmask = false(clustersize,clustersize,nclusters);
    
    for cl = 1 : nclusters
        y = ceil(cl/nclx);
        x = cl - (y-1)*nclx;
        j = (x-1)*clustersize + 1 : x*clustersize;
        i = (y-1)*clustersize + 1 : y*clustersize;
        cllandcover(:,:,:,cl) = landcover(i,j,:);
        clecoregion(:,:,cl) = ecoregion(i,j);
        clpixarea(:,:,cl) = pixelarea(i,j);
        cllandmask(:,:,cl) = landmask(i,j);
    end
    
    parfor cl = 1 : nclusters
        cllc = cllandcover(:,:,:,cl);
        cleco = clecoregion(:,:,cl);
        clpixarea = clecoregion(:,:,cl);
        clprevlc = clprevalentlc(:,:,cl);
        mask = cllandmask(:,:,cl);
        
        for tl = 1 : ntlpcl
            y = ceil(tl/ntlx);
            x = tl - (y-1)*ntlx;
            j = (x-1)*tilesize + 1 : x*tilesize;
            i = (y-1)*tilesize + 1 : y*tilesize;
                        
            smmask = mask(i,j);
            
            if sum(smmask,'all') > 0
                smmodis = cllc(i,j,:);
                smprev = clprevlc(i,j);
                smpix = clpixarea(i,j);
                smeco = cleco(i,j);
                ecolist = setdiff(unique(smeco),econoval);
                
                for ee = 1 : numel(ecolist)
%                     flag = true;
                    ii = smeco == ecolist(ee);
                    ecoarea = sum(smpix(ii),'omitnan');
%                     if ecoarea/100 < minecopixels, flag = false; end
                    ta = zeros(numel(landcoverindex),1); ct = 1;
                    while max(ta) == 0 && ct <= numel(modisvarnames)                        
                        smlandcover = smmodis(:,:,ct);
                        for ll = 1 : numel(landcoverindex)
                            jj = smlandcover == landcoverindex(ll);
                            kk = ii + jj == 2;
                            ta(ll) = sum(smpix(kk),'all','omitnan');
                        end
%                         if flag == true
                            validta = (ta ./ repmat(ecoarea,[numel(landcoverindex),1]))*100 >= ...
                                prevalentlandcoverthreshold;
%                         else
%                             validta = ta >= minecopixels;
%                         end
                        ta(logical(1-validta)) = 0;
                        ct = ct + 1;
                    end
                    if max(ta) > 0
                        m = find(ta==max(ta));
                        if numel(m) == 1
                            smprev(ii) = landcoverindex(m);
                        elseif ismember(m,mostlikelypereco(ecolist(ee))) %#ok<PFBNS>
                            smprev(ii) = mostlikelypereco(ecolist(ee));
                        end
                    end
                end
                
                clprevlc(i,j) = smprev;
            end            
        end
        
        clprevalentlc(:,:,cl) = clprevlc;
    end
    for cl = 1 : nclusters
        y = ceil(cl/nclx);
        x = cl - (y-1)*nclx;
        j = (x-1)*clustersize + 1 : x*clustersize;
        i = (y-1)*clustersize + 1 : y*clustersize;
        prevalentlc(i,j) = clprevalentlc(:,:,cl);
    end

    
else
    
    for tl = 1 : ntiles
        y = ceil(tl/nblx);
        x = tl - (y-1)*nblx;
        j = (x-1)*tilesize + 1 : x*tilesize;
        i = (y-1)*tilesize + 1 : y*tilesize;

        smmask = landmask(i,j);

        if sum(smmask,'all') > 0
            smmodis = landcover(i,j,:);
            smprev = prevalentlc(i,j);
            smpix = pixelarea(i,j);
            smeco = ecoregion(i,j);
            ecolist = setdiff(unique(smeco),econoval);
            
            for ee = 1 : numel(ecolist)
%                 flag = true;
                ii = smeco == ecolist(ee);
                ecoarea = sum(smpix(ii),'omitnan');
%                 if ecoarea/100 < minecopixels, flag = false; end
                ta = zeros(numel(landcoverindex),1); ct = 1;
                while max(ta) == 0 && ct <= numel(modisvarnames)
                    smlandcover = smmodis(:,:,ct);
                    for ll = 1 : numel(landcoverindex)
                        jj = smlandcover == landcoverindex(ll);
                        kk = ii + jj == 2;
                        ta(ll) = sum(smpix(kk),'all','omitnan');
                    end
%                     if flag == true
                        validta = (ta ./ repmat(ecoarea,[numel(landcoverindex),1]))*100 >= ...
                            prevalentlandcoverthreshold;
%                     else
%                         validta = ta >= minecopixels;
%                     end
                    ta(logical(1-validta)) = 0;
                    ct = ct + 1;
                end
                if max(ta) > 0
                    m = find(ta==max(ta));
                    if numel(m) == 1
                        smprev(ii) = landcoverindex(m);
                    elseif ismember(m,mostlikelypereco(ecolist(ee)))
                        smprev(ii) = mostlikelypereco(ecolist(ee));
                    end
                end
            end
                        
            prevalentlc(i,j) = smprev;
        end
    end

end

prevalentlc(landmask==false) = u8noval;