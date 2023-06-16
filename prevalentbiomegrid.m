function [biome,biobl,rlst,clst] = prevalentbiomegrid(resolution,landproportion,region,...
    climate,area,RegionalMostLikelyBiome,landmask,biomeindex,latlonscale,novalue)

% FUNCTION prevalentbiomegrid: find prevalent biome type in a grid cell
%   For each grid cell, this function determines the prevalent biome-type in all land pixels within
%   the grid cell that belong to the same region and climate zone. The grid cell is of a different
%   resolution ("resolution") than the original land cover map of land cover proportion (which is in
%   the "latlonscale" resolution). Output map is in the original resolution "latlonscale".
% Input variables are as follows:
%   1.  resolution                  Resolution of the "grid" map or grid cell resolution
%   2.  landproportion              [nlat,nlon,nlc] land cover proportion array of all nlc land
%                                       cover of interest given in "native" resolution "latlonscale"
%   3.  region                      [nlat,nlon] global regions array
%   4.  climate                     [nlat,nlon] global climate zones array
%   5.  area                        [nlat,nlon] pixel area array
%   6.  RegionalMostLikelyBiome     [nreg,nclim] array of most prevalent biome (in a given
%                                       biome-of-interest list) for each region (lines) and climate
%                                       (columns) combination
%   7.  landmask                    [nlat,nlon] global array of land pixels with values of
%                                       0 - water
%                                       1 - land
%   8.  biomindex                   list of biome of interest indexes
%   9.  latlonscale                 resolution of the original maps. Default is 0.05 decimal degrees
%   10. novalue                     no-data value for input/output maps. Default is 255


if nargin < 10
    novalue = 2^8-1;
end
if nargin < 9
    latlonscale = 0.05;
end

[blocksize,~,nbt] = size(landproportion);

if nargin < 8
    biomeindex = 1:nbt;
end
if nargin < 7
    landmask = true(blocksize,blocksize);
    landmask(region==novalue) = false;
end

tilesize = resolution / latlonscale;
nblx = blocksize / tilesize;
nbly = blocksize / tilesize;
ntiles = nblx * nbly;

biome = ones(blocksize,blocksize,'uint8') .* novalue;
biome(landmask) = 0;

biobl = [];
rlst = [];
clst = [];

for tl = 1 : ntiles
    y = ceil(tl/nblx);
    x = tl - (y-1)*nblx;
    j = (x-1)*tilesize + 1 : x*tilesize;
    i = (y-1)*tilesize + 1 : y*tilesize;
    
    tilelandmask = landmask(i,j);
    tilepf = ones(tilesize,tilesize,'uint8') * novalue;
    tilepf(tilelandmask) = 0;
    pbio = landproportion(i,j,:);
    
    if sum(pbio,'all') > 0
        tileclim = climate(i,j);
        tilereg = region(i,j);
        tilearea = area(i,j);
        
        rlst = setdiff(unique(tilereg),novalue);
        clst = setdiff(unique(tileclim),novalue);
        
        if ntiles == 1
            biobl = zeros(numel(rlst),numel(clst));
        end
        
        for rr = 1 : numel(rlst)
            for cc = 1 : numel(clst)
                thisreg = rlst(rr);
                thisclim = clst(cc);
                k = ((tilereg == thisreg) + (tileclim == thisclim)) == 2;
                mlf = RegionalMostLikelyBiome(thisreg,thisclim);
                
                A = zeros(nbt,1);
                for bb = 1 : nbt
                    bp = pbio(:,:,bb);
                    A(bb) = sum(double(bp(k)) .* tilearea(k));
                end
                
                if sum(A) == 0
                    val = 0;
                elseif sum(A == max(A)) == 1
                    val = find(A==max(A),1);
                else
                    f2 = find(A==max(A));
                    if ismember(mlf,f2)
                        val = mlf;
                    else
                        val = f2(random('Discrete Uniform',numel(f2)));                       
                    end
                end
                if val > 0
                    tilepf(k) = biomeindex(val);
                    if numel(biobl) > 0
                        biobl(rr,cc) = biomeindex(val); %#ok<AGROW>
                    end
                end
            end
        end
        biome(i,j) = tilepf;
    end    
    
end
