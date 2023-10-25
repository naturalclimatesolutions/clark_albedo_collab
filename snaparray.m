function [newmaparray] = snaparray(refmap,refoldmap,oldmap)

% Get map reference type and no data values
cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double"];
mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan];

if isa(refmap,'map.io.RasterInfo'), refmap = refmap.RasterReference; end
if isa(refoldmap,'map.io.RasterInfo')
    noval = refoldmap.MissingDataIndicator;
    if numel(noval) == 0
        fcl = refoldmap.NativeFormat;
        ii = ismember(cllist,fcl);
        noval = mislist(ii);
    end
    refoldmap = refoldmap.RasterReference; 
else
    fcl = class(oldmap);
    ii = ismember(cllist,fcl);
    noval = mislist(ii);
end

% Get pixel size and map sizes
mapsize = refmap.RasterSize;
oldmapsize = refoldmap.RasterSize;
if strcmp(refmap.CoordinateSystemType,'geographic')
    pixelsize = refmap.CellExtentInLatitude;
    oldpixelsize = refoldmap.CellExtentInLatitude;
    reflonlimits = refmap.LongitudeLimits;
    oldlonlimits = refoldmap.LongitudeLimits;
    reflatlimits = refmap.LatitudeLimits;
    oldlatlimits = refoldmap.LatitudeLimits;
elseif strcmp(refmap.CoordinateSystemType,'planar')
    pixelsize = refmap.CellExtentInWorldX;
    oldpixelsize = refoldmap.CellExtentInWorldX;
    reflonlimits = refmap.XWorldLimits;
    oldlonlimits = refoldmap.XWorldLimits;
    reflatlimits = refmap.YWorldLimits;
    oldlatlimits = refoldmap.YWorldLimits;
else
    error("code for your CoordinateSystemType")
end

if abs(oldpixelsize - pixelsize) > 10^-6
    error("the map to be snapped has the wrong pixel size!!");
end

% Get no data value (when adding grid cells)

% Get number of pixel differences
ftdx = round((reflonlimits - oldlonlimits) ./ pixelsize);
ftdy = round((reflatlimits - oldlatlimits) ./ pixelsize);

% Because I have to round, I might be off by one pixel or so, correct if it's the case
nsx = oldmapsize(2) - ftdx(1) + ftdx(2);
dum = abs(ftdx - ((reflonlimits - oldlonlimits) ./ pixelsize));
% case when I take out one pixel too many in x
if nsx < mapsize(2) && ftdx(1) > 0
    if dum(1) > dum(2), ftdx(1) = ftdx(1) - 1; else, ftdx(2) = ftdx(2) + 1; end
    % case when I don't add enough pixels in x
elseif nsx < mapsize(2) && ftdx(1) < 0
    if dum(1) > dum(2), ftdx(2) = ftdx(2) + 1; else, ftdx(1) = ftdx(1) - 1; end
    % Case when I add one too many pixels in x
elseif nsx > mapsize(2) && ftdx(1) < 0
    if dum(1) > dum(2), ftdx(1) = ftdx(1) + 1; else, ftdx(2) = ftdx(2) - 1; end
    % Case when I don't take out enough pixels in x
elseif nsx > mapsize(2) && ftdx(1) > 0
    if dum(1) > dum(2), ftdx(2) = ftdx(2) - 1; else, ftdx(1) = ftdx(1) + 1; end
end
nsx = oldmapsize(2) - ftdx(1) + ftdx(2);
if nsx ~= mapsize(2), error("I did not corect my ftdx array properly"); end

nsy = oldmapsize(1) - ftdy(1) + ftdy(2);
dum = abs(ftdy - ((reflatlimits - oldlatlimits) ./ pixelsize));
% case when I take out one pixel too many in x
if nsy < mapsize(1) && ftdy(1) > 0
    if dum(1) > dum(2), ftdy(1) = ftdy(1) - 1; else, ftdy(2) = ftdy(2) + 1; end
    % case when I don't add enough pixels in x
elseif nsy < mapsize(1) && ftdy(1) < 0
    if dum(1) > dum(2), ftdy(2) = ftdy(2) + 1; else, ftdy(1) = ftdy(1) - 1; end
    % Case when I add one too many pixels in x
elseif nsy > mapsize(1) && ftdy(1) < 0
    if dum(1) > dum(2), ftdy(1) = ftdy(1) + 1; else, ftdy(2) = ftdy(2) - 1; end
    % Case when I don't take out enough pixels in x
elseif nsy > mapsize(1) && ftdy(1) > 0
    if dum(1) > dum(2), ftdy(2) = ftdy(2) - 1; else, ftdy(1) = ftdy(1) + 1; end
end
nsy = oldmapsize(1) - ftdy(1) + ftdy(2);
if nsy ~= mapsize(1), error("I did not corect my ftdy array properly"); end


    

% Modify map in the x-direction (West-East)
dum = oldmap;
if ftdx(1) < 0
    adx1 = ones(oldmapsize(1),abs(ftdx(1))) .* noval;
    dum = cat(2,adx1,dum);
elseif ftdx(1) > 0
%     warning(strcat("reference map starts further east than map that needs adjusting, ",...
%         "see if it is a problem!"))
    dum = dum(:,1+ftdx(1):oldmapsize(2));
end
oldmapsize(2) = oldmapsize(2) - ftdx(1);
if ftdx(2) > 0
    adx2 = ones(oldmapsize(1),abs(ftdx(2))) .* noval;
    dum = cat(2,dum,adx2);
elseif ftdx(2) < 0
%     warning(strcat("reference map ends further west than map that needs adjusting, ",...
%         "see if it is a problem!"))   
    dum = dum(:,1:oldmapsize(2)+ftdx(2));
end
oldmapsize(2) = oldmapsize(2) + ftdx(2);

% Modify map in the y-direction. Given that smaller latitudes are further south, the algorithm is
% reversed
if ftdy(2) > 0
    ady2 = ones(ftdy(2),oldmapsize(2)) .* noval;
    dum = cat(1,ady2,dum);
elseif ftdy(2) < 0
%     warning(strcat("reference map ends further north than map that needs adjusting, ",...
%         "see if it is a problem!"))   
    dum = dum(1-ftdy(2):oldmapsize(1),:);
end
oldmapsize(1) = oldmapsize(1) + ftdy(2);
if ftdy(1) < 0
    ady1 = ones(abs(ftdy(1)),oldmapsize(2)) .* noval;
    dum = cat(1,dum,ady1);
elseif ftdy(1) > 0
%     warning(strcat("reference map starts further south than map that needs adjusting, ",...
%         "see if it is a problem!"))   
    dum = dum(1:oldmapsize(1)-ftdy(1),:);
end
oldmapsize(1) = oldmapsize(1) - ftdy(1);

% Save new map
newmaparray = dum;

% Make sure I did it right ...
if isempty(find(oldmapsize - size(dum),1)) == 0
    warning('I did not add my sizes correctly in my map correction')
    return
end
if isempty(find(oldmapsize - mapsize,1)) == 0
    warning('I did not resize my forest map correctly')
    return
end
if sum(size(newmaparray) - mapsize) ~= 0
    warning('I did not resize my forest map correctly')
    return
end
    

