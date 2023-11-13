function [lowresmap,highresmap] = printmapprep(input,fileroot,datablock,inputresolution,...
    inputdataformat,inputsize,blocksize,lowerresolution,aggregationmethod,missinglandblock,...
    indexfileroot,mask,flag)

% PRINTMAPPREP
%   This function is used to create two arrays, both representing the entire geographical map of a
%   given variable. The "highresmap" is the high resolution array used to export the map as GeoTIFF
%   format. The "lowresmap" is the decreased resolution used to plot the variable in Matlab as the
%   high resolution is too big for Matlab graphics to process.
%
% Input variables are
%   - input                 % name of input variable, Variable can either be a single variable
%                               (f.eg. albedo) or a structure array. The structure array is
%                               used if the desired variable is the difference or the ratio between
%                               two variables (f.eg. difference between NCI when using Walker or 
%                               ESA-truncated Walker). Structure needs to be passed as follows                               
%                                   input.base - base variable (f.eg. NCI-Walker)
%                                   input.modified  - compared variable (f.eg. NCI-truncated ESA)
%                                   input.method -  comparison method (choices are difference,
%                                                   percent and proportion)
%   - fileroot              % root name of subdomain files where the input variable is stored
%   - datablock             % logical array describing presence or absence of block file
%   - inputresolution       % resolution of input variable (high resolution)
%   - inputdataformat       % data format of input variable
%   - inputsize             % size of high resolution map array
%   - blocksize             % size of variable as saved in subdomain files
%   - lowerresolution       % resolution of coarse array used for plotting
%   - aggregationmethod     % method used to aggregate data in lower resolution, choices are
%                               mean or median for continuous data and mode for categorical data
%   - missinglandblock      % subdomain block number that contains land area but that was not saved
%                               as a file because there was no input data in it. F.eg. block #51
%                               (in a 10x10 subdivision) is a block in Greenland that had no data in
%                               the Walker dataset and thus is missing from the file list
%   - indexfileroot         % subdomain file root name where x-y indeces are stored. By default it's
%                               the same as fileroot
%   - mask                  % name of a specific mask to be applied to the input data. This is used
%                               when for opportunity maps, where we only select specific locations
%                               on the globe. By default, no mask is applied.
%   - flag                  % (optional) flag is used to decide if 0 values are to be considered as
%                               nodata value (flag=true) or kept as zero values (flag=false).
%                               Default is false.
%
% Created by Natalia Hasler, Clark University - April 2023


if nargin < 13
    flag = false;
end
if nargin < 12
    mask = 1;
    maskdata = true([blocksize,blocksize]);
end
if nargin < 11
    indexfileroot = fileroot;
end
if nargin < 10
    missinglandblock = [];
end
if nargin < 9
    aggregationmethod = "mean";
end
if isa(input,'struct')
    inpar = string(fields(input));
    if sum(ismember(inpar,["base","modified","method"])) ~= 3
        error("The input structure passed into the function does not have the expected fields")
    end
    varbase = input.base;
    varmod = input.modified;
    vmethod = input.method;
    invars = [varmod,varbase];
else
    invars = input;
    inpar = "base";
end

cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double","logical"];
mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan,false];
nodatavalue = mislist(cllist==inputdataformat);
ratio = lowerresolution / inputresolution;
lowerressize = inputsize ./ ratio;
nblx = inputsize(2) ./ blocksize;
lowresblocksize = round(blocksize ./ ratio);    % rounding is to avoid warnings

highresmap = ones(inputsize,inputdataformat) .* nodatavalue;
lowresmap = ones(lowerressize,inputdataformat) .* nodatavalue;

if inputresolution == 0.005
    indexvars = ["x","y","smregi","smregj"];
else
    indexvars = ["x","y","regi","regj"];
end

if strcmp(mask,"walkeroppcat")
    Wsc = [2,5];
    warning(strcat("For ease of use I hardcoded the Walker Opportunity classes as being #2 & #5,",...
        " modify code if that is not the case anymore!"))
end


for bb = 1 : numel(datablock)
    if datablock(bb) == false
        if bb==missinglandblock && (~isempty(missinglandblock) && exist('misblockvalhr','var'))
            y = ceil(bb/nblx);
            x = bb - (y-1)*nblx;
            smregi = (y-1)*blocksize + 1 : y*blocksize;
            smregj = (x-1)*blocksize + 1 : x*blocksize;
            highresmap(smregi,smregj) = misblockvalhr;
            regj = (x-1)*lowresblocksize + 1 : x * lowresblocksize;
            regi = (y-1)*lowresblocksize + 1 : y * lowresblocksize;
            lowresmap(regi,regj) = misblockvallr;
            clear smregi smregj regi regj x y 
        end
        continue
    end
    subdatafname = strcat(fileroot,num2str(bb),".mat");
    indfname = strcat(indexfileroot,num2str(bb),".mat");
    load(indfname,indexvars{:})
    if numel(inpar) == 1
        load(subdatafname,invars)
        if ~exist(invars,'var'), continue; end
        eval(strcat("dum = ",invars,";"));
        if numel(dum) == 0, continue; end
    else
        load(subdatafname,invars{:})
        eval(strcat("dum = ",invars(1),";"));
        if numel(dum) == 0, continue; end
    end
    
    if ~isnumeric(mask)
        load(indfname,mask)
        switch mask
            case "walkeroppcat"
                maskdata = ismember(walkeroppcat,Wsc);
            case "griscom"
                maskdata = griscom == 1;
            otherwise
                eval(strcat("maskdata =",mask,";"))
        end
    end
    
    if exist("regi",'var')
        smregi = regi;
        smregj = regj;
        clear regi regj
    end
    
    if numel(inpar) == 1
        eval(strcat("data = ",input,";"));
    else
        eval(strcat("vbase = ",varbase,";"));
        eval(strcat("vmod = ",varmod,";"));
        switch vmethod
            case "difference"
                data = vmod - vbase;
            case "percent"
                data = (vmod - vbase) ./ vbase .*100;
            case "proportion"
                data = (vmod - vbase) ./ vbase;
            otherwise
                error("Do not know this method")
        end
    end
    
    dataclass = class(data);
    elimdata = logical(1-maskdata);
    data(elimdata) = mislist(cllist==dataclass); 
    
    if ismember(dataclass,["single","double"]) && flag == true
        data(data==0) = nan;  
    end
    
    if ~strcmp(dataclass,inputdataformat)
        if ismember(dataclass,["single","double"]) % Allowing a change in format for size issues
            data(isnan(data)) = nodatavalue;  
            eval(strcat("hrdata =",inputdataformat,"(data);"))
        else
            error("Inputdataformat is different from data format of input values")
        end
    else
        hrdata = data;
    end

    highresmap(smregi,smregj) = hrdata;
    
    if ~ismember(dataclass,["single","double","logical"])
        data = single(data);
        data(data==nodatavalue) = NaN;
    end
    
    regj = (x-1)*lowresblocksize + 1 : x * lowresblocksize;
    regi = (y-1)*lowresblocksize + 1 : y * lowresblocksize;
    smallmap = nan(lowresblocksize,lowresblocksize);
    
    for i = 1 : lowresblocksize
        ii = (i-1)*ratio + 1 : i * ratio;
        for j = 1 : lowresblocksize
            jj = (j-1)*ratio + 1 : j * ratio;
            switch aggregationmethod
                case "mean"
                    val = mean(data(ii,jj),'all','omitnan');
                case "mode"
                    val = mode(data(ii,jj),'all');
                case "median"
                    val = median(data(ii,jj),'all','omitnan');
            end
            if isnan(val)
                smallmap(i,j) = nodatavalue;
            else
                smallmap(i,j) = val;
            end
        end
    end
    
    lowresmap(regi,regj) = smallmap;

    if bb == missinglandblock-1
        misblockvalhr = data;
        misblockvallr = smallmap;
    end

    if numel(inpar) == 1
        clear(invars,input)
    else
        clear(invars{:})
    end
    clear(indexvars{:})
    clear regi regj smallmap i ii j jj xsize ysize data hrdata vbase vmod
    if ceil(bb/50) == bb/50
    strcat("Done with map preparation for ",invars(1)," in block #",num2str(bb))
    end
        
end
