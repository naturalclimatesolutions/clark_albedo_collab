function blc = bestlandcover(latindex,lonindex,landcoverproportion,landcoverlist,...
    climate,region,area)

% FUNCTION bestlandcover: defines land cover (from a defined list) with most area in a pixel
%   This function assigns the land cover class with the most area to a given pixel. Land cover class
%   is chosen amoung a pre-defined subset of land cover classes given in "landcoverlist". When
%   several land cover classes have the same area, class with most area in the 9x9 cluster centered
%   aroung pixel is chosen. Only pixels with same region and climate indices as the pixel of
%   interest are considered in the cluster.
% Input variables are as follows:
%   1.  latindex                [npix] latitude index of one or several pixels (npix)
%   2.  lonindex                [npix] longitude index of the same pixel(s)
%   3.  landcoverproportion     [nlat,nlon,nlc] land cover proportions of all nlc number of land
%                                   cover of interest, with nlat and nlon the number of indices in
%                                   latitude/longitude to cover the globe
%   4. landcoverlist            list of class numbers of all land cover of interest
%   5. climate                  [nlat,nlon] global climate zones array
%   6. region                   [nlat,nlon] global region array
%   7. area                     [nlat,nlon] area of each pixel globally


if numel(latindex) ~= numel(lonindex)
    error("please provide pixel coordinates in [latindex,lonindex] !")
end

nlon = 7200;
[sizey,sizex] = size(climate);
blc = zeros(numel(latindex),1);

for px = 1 : numel(latindex)
    a = latindex(px);
    b = lonindex(px);
    prop = landcoverproportion(a,b,:);
    clim = climate(a,b);
    reg = region(a,b);
    if sum(prop>0)>0
        if sum(prop == max(prop)) == 1
            blc(px) = landcoverlist(find(prop==max(prop),1));
        else
            smi = a-1:a+1; smi(smi<1) = NaN; smi(smi>sizey) = NaN; smi = smi(isfinite(smi));
            smj = b-1:b+1;
            if sizex == nlon
                smj(smj<1) = smj(smj<1)+nlon; smj(smj>sizex) = smj(smj>sizex)-nlon;
            else
                smj(smj<1) = NaN; smj(smj>sizex) = NaN; smj = smj(isfinite(smj));
            end
            sma = area(smi,smj);
            smp = landcoverproportion(smi,smj,:);
            scl = climate(smi,smj) == clim;
            sre = region(smi,smj) == reg;
            vpx = scl + sre == 2;
            p2 = zeros(numel(landcoverlist),1);
            for lc = 1 : numel(landcoverlist)
                pp = smp(:,:,lc);
                p2(lc) = sum(double(pp(vpx)) .* sma(vpx));
            end
            k = landcoverlist(p2 == max(p2));
            if numel(k) == 1
                blc(px) = k;
            else
                blc(px) = k(random('Discrete Uniform',numel(k)));
            end
        end
    end
end