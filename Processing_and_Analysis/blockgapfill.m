function outputvar = blockgapfill(inputvar,mask,tilesize,nodatavalues,maxstep,method)

% FUNCTION Blockgapfill: Geographical "categorical" gap filling over a geographical subset (block)
%   This function is used to gap-fill a geo-referenced categorical variable with nearest neighbors.
%   It is really only intended for variable such as climate zones or ecoregions, which can have a
%   few missing pixels here and there due to discrepancies in the definition of land masks. It is 
%   similar to the geogapfill function, but used on a geographical subset only, tipically a tile of
%   10x10 deg lat-lon or 20x20 latlon. Because it uses the prevalent value (mode), it is not meant
%   for continuous variables datasets.
%
% Input variables are:
%   inputvar                [ny,nx] array to be gap-filled, also allowed [ny,nx,nm]
%   mask                    [ny,nx] logical array of land data points of interest
%   tilesize                nx/ny dimension (here it's a square)
%   nodatavalues            list of nodata values in the inputvar array (first one is the one that
%                               will be prioritized as nodata value in output)
%   maxstep                 [optional] maximum number of pixels in the "search zone" for
%                               gap-filling. Default is 50
%
% Output variable is same as inputvar array except for gap-filled pixels
%
% Created by Natalia Hasler, Clark University - January 2023


% Make sure there is data available
vals = unique(inputvar);
outputvar = inputvar;

if numel(vals) == 1 && ismember(vals,nodatavalues)
    warning("There is no availaible data in the input file")
    return
end

% Define variables
if nargin < 6
    method = "mode";
end
if nargin < 5
    maxstep = 50;
end
if nargin < 4
    cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double","logical"];
    mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan,false];
    dataclass = class(inputvar);
    nodatavalues = mislist(strcmp(cllist,dataclass));
end
if nargin < 3
    [tilesize,~] = size(inputvar);
end

% Prepare temporary variables to identify missing data
tmpvar = double(inputvar);
for k = 1 : numel(nodatavalues)
    ndv = nodatavalues(k);
    tmpvar(tmpvar==ndv) = nan;
end
tmpvar2 = tmpvar;
tmpvar(mask==0) = -1;

% Loop through missing values
pts = find(isnan(tmpvar));
[a,b,m] = ind2sub(size(tmpvar),pts);
for pt = 1 : numel(a)
    ii = a(pt);
    jj = b(pt);
    mm = m(pt);
    val = tmpvar(ii,jj,mm);
    st = 0;
    while isnan(val) && st <= maxstep % When no values are found within +-'maxstep', ignore them
        st = st + 1;
        smi = ii-st:ii+st; smi(smi<1) = 1; smi(smi>tilesize) = tilesize;
        smj = jj-st:jj+st; smj(smj<1) = 1; smj(smj>tilesize) = tilesize;
        tmpvar(smi,smj,mm) = tmpvar2(smi,smj,mm);
        switch method
            case "mode"
                val = mode(tmpvar(smi,smj,mm),'all');
            case "mean"
                val = mean(tmpvar(smi,smj,mm),'all','omitnan');
        end
    end
    if ~isnan(val)
        tmpvar(ii,jj) = val; tmpvar2(ii,jj) = val;
    end
%     if floor(pt/1000) == pt/1000
%         strcat("Done with pt #",num2str(pt))
%     end
end

% save gap-filled values
tmpvar2(isnan(tmpvar2)) = nodatavalues(1);
outputvar(mask) = tmpvar2(mask);
