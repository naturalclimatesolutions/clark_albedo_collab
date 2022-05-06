
function outputvar = geogapfill(latitude,longitude,inputvar,varname,nanval,globalpixelarea,...
    latitudinalprevalence,region)

% FUNCTION Geographical gap fill
%   This function is used to gap-fill a geo-referenced variable with nearest neighbors. It is really
%   only intended for variable such as climate zones or regions, which can have a few pixels
%   missing here and there due to discrepancies in the definition of land masks. This occurs on
%   coastal regions or on very small islands ... It also occurs in Antarctica where MODIS landmask
%   differs from other landmasks. I use first the nearest neighbors within 50 pixels (in each
%   direction). Outside that distance, I define crude boundaries of different regions to assign
%   a value to small islands. As for climate zones, I will use the most prevalent in that latitude
%   as long as that climate exists in the region, otherwise I'll again find most likely given
%   latitude and region. For gap fillings with different constraints, edit the code.
%
% Input variables are:
%   latitude                latitude index in the inputvar array
%   longitude               longitude index in the inputvar array
%   inputvar                [nlat,nlon] array to be gap-filled
%   varname                 string array with the name list of the different values in inputvar, 
%                               those will be the region names or the climate names
%   nanval                  value of nodata in the inputvar array
%   globalpixelarea         [optional] area for each lat/lon pixel on the globe
%   latitudinalprevalence   [optional] most common value of inputvar in any given latitude
%   region                  [optional] when inputvar is not "regions", defines region value for each
%                               land pixel. This is used to calculate so that gap-filled value is
%                               consistent with values in region.
%
% Output variable is same as inputvar array except for the additional [latitude,longitude] pixel
%   that has been gap-filled

[nlat,nlon] = size(inputvar);
latlonscale = 180 / nlat;
latvals = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lonvals = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;

lat = latvals(latitude);
lon = lonvals(longitude);

if nargin < 6
    earthellipsoid = referenceSphere('earth','m');
    pix = zeros(nlat,1);
    lat1 = 90;
    lon1 = -180;
    lon2 = lon1 + latlonscale;
    for i = 1 : nlat
        lat2 = lat1 - latlonscale;
        pix(i) = areaquad(lat1,lon1,lat2,lon2,earthellipsoid);
        lat1 = lat2;
    end
    globalpixelarea = repmat(pix,[1 nlon]);
end
if nargin < 7
    latitudinalprevalence = [];
else
    if isa(latitudinalprevalence,'uint8')
        latitudinalprevalence = double(latitudinalprevalence);
        latitudinalprevalence(latitudinalprevalence==2^8-1) = NaN;
    else
        error("Code here for a different format of latitudinalprevalence")
    end
end
if nargin < 8
    region = [];
elseif isnan(region(1,1))
    reglist = unique(region(isfinite(region)));
else
    reglist = unique(region);
    reglist = setdiff(reglist,region(1,1));
end
if sum(strcmpi("South America",varname)) == 1
    regionnames = varname;
end

Polynesia = [33,-180;26,-157;-28,-109;-47,-180;0,-180];
Micronesia = [0,130;23,134;25,165;-6,180];
Melanesia = [0,130;-6,180;-34,180;-30,155];
CA = [7,-78;7,-110;30,-130;32.5,-117;25,-97];
SA = [7,-78;7,-100;-55,-90;-55,-45;-5,-15;7,-55];
NA = [82,-8;70,-19;60,-40;30,-40;30,-80;60,-80];
NE = [82,-8;70,-19;60,-30;50,-5;70,30;90,30];
MA = [-17,12;-17,-20;4,-20;4,12];
SAf = [-17,10;-17,-20;-55,-20;-55,45;-30,30];
ERus = [69,31;85,35;85,73;69,64];
EA = [-27,32;-45,85;12,64;12,45];
SAs=[21,92;6.5,95;-10,78;25,60];
SEA = [-8.5,137;-20,90;6,94.5;16,95;21,121;-1,138];
EAs = [18,110;22,140;43,160;45,142];
Aus = [-10,142;-40,179;-55,160;-50,100;-20,107];

% Russia = [59,-180;66,-169;90,-180];
% Caribbean = [30,-80;18,-55;9,-60;8,-82;18,-94;26,-97]; %%

outputvar = inputvar;
step = 0; cond = 0;
val = inputvar(latitude,longitude);

while cond == 0
    step = step + 1;
    i = latitude-step:latitude+step; i(i<1) = NaN; i(i>nlat) = NaN; i = i(isfinite(i));
    j = longitude-step:longitude+step; j(j<1) = j(j<1)+nlon; j(j>nlon) = j(j>nlon)-nlon;
    zonevar = double(inputvar(i,j));
    if isfinite(nanval)
        zonevar(zonevar==nanval) = NaN;
    end
    if sum(isfinite(zonevar),'all')>0
        s = min(zonevar,[],'all','omitnan'):max(zonevar,[],'all','omitnan');
        n = histcounts(zonevar);
        if sum(n == max(n)) == 1
            val = s(n==max(n));
            cond = 1;
        else
            A = zeros(length(s),1);
            for k = 1 : length(s)
                pp = zonevar == s(k);
                A(k) = sum(globalpixelarea(pp));
            end
            kk = find(A == max(A));
            if numel(kk) == 1
                val = s(kk);
                cond = 1;
            end
        end
    elseif latvals(latitude) <= -60
        if sum(ismember(varname,"Antarctica")) == 1
            val = find(ismember(varname,"Antarctica"),1);
        elseif sum(ismember(varname,"SNO")) == 1
            val = find(ismember(varname,"SNO"),1);
        elseif sum(ismember(varname,"South America")) == 1
            val = numel(varname) + 1;
            warning(strcat("region/continent did not include Antarctica, so I added value ",...
                num2str(numel(varname)+1)," as the value for Antarctica"))
        else
            error("I have missing values in Antartica, and don't know how to fill them!")
        end
        cond = 1;
    elseif step > 50
        if numel(latitudinalprevalence) > 0
            [nllp,nreg] = size(latitudinalprevalence);
            chunk = nlat/nllp;
            if numel(region) > 0, reg = region(latitude,longitude); else, reg = NaN; end
            lati = ceil(latitude/chunk);
            if isnan(reg) || nreg == 1
                val = latitudinalprevalence(lati,nreg);
            elseif nreg == numel(reglist) + 1
                val = latitudinalprevalence(lati,reg);
            else
                error("The 'latitudinalprevalence' array does not have the expected size!")
            end
            if isnan(val)
                list = unique(latitudinalprevalence(:,reg));
                list = list(isfinite(list));
                if ismember(latitudinalprevalence(lati,nreg),list)
                    val = latitudinalprevalence(lati,nreg);
                else
                    cond2 = 0; step2 = 0;
                    if isnan(reg), b = nreg; else, b = [reg,nreg]; end
                    while cond2 == 0
                        step2 = step2 + 1;
                        a = lati-step2:lati+step2; a = a(a>0); a = a(a<=nllp);
                        dum = latitudinalprevalence(a,b);
                        dum2 = intersect(dum,list);
                        if numel(dum2) == 1
                            val = dum2;
                            cond2 = 1;
                        elseif numel(dum2) > 1
                            s = dum2(1):dum2(numel(dum2));
                            n = histcounts(dum);
                            if sum(n == max(n)) == 1
                                val = s(n==max(n));
                                cond2 = 1;
                            else
                                dum3 = find(n==max(n));
                                val = dum3(random('Discrete Uniform',numel(dum3)));
                                cond2 = 1;
                            end
                        end
                    end
                end
            end
        else
            in = inpolygon(lat,lon,Polynesia(:,1),Polynesia(:,2));
            if in
                val = find(strcmpi("Polynesia",regionnames));
            else
                in = inpolygon(lat,lon,Micronesia(:,1),Micronesia(:,2));
                if in
                    val = find(strcmpi("Micronesia",regionnames));
                else
                    in = inpolygon(lat,lon,Melanesia(:,1),Melanesia(:,2));
                    if in
                        val = find(strcmpi("Melanesia",regionnames));
                    else
                        in = inpolygon(lat,lon,CA(:,1),CA(:,2));
                        if in
                            val = find(strcmpi("Central America",regionnames));
                        else
                            in = inpolygon(lat,lon,SA(:,1),SA(:,2));
                            if in
                                val = find(strcmpi("South America",regionnames));
                            else
                                in = inpolygon(lat,lon,NA(:,1),NA(:,2));
                                if in
                                    val = find(strcmpi("Northern America",regionnames));
                                else
                                    in = inpolygon(lat,lon,NE(:,1),NE(:,2));
                                    if in
                                        val = find(strcmpi("Northern Europe",regionnames));
                                    else
                                        in = inpolygon(lat,lon,MA(:,1),MA(:,2));
                                        if in
                                            val = find(strcmpi("Middle Africa",regionnames));
                                        else
                                            in = inpolygon(lat,lon,SAf(:,1),SAf(:,2));
                                            if in
                                                val = find(strcmpi("Southern Africa",regionnames));
                                            else
                                                in = inpolygon(lat,lon,ERus(:,1),ERus(:,2));
                                                if in
                                                    val = find(strcmpi("European Russia",regionnames));
                                                else
                                                    in = inpolygon(lat,lon,EA(:,1),EA(:,2));
                                                    if in
                                                        val = find(strcmpi("Eastern Africa",regionnames));
                                                    else
                                                        in = inpolygon(lat,lon,SAs(:,1),SAs(:,2));
                                                        if in
                                                            val = find(strcmpi("Southern Asia",regionnames));
                                                        else
                                                            in = inpolygon(lat,lon,SEA(:,1),SEA(:,2));
                                                            if in
                                                                val = find(strcmpi("Southeastern Asia",regionnames));
                                                            else
                                                                in = inpolygon(lat,lon,EAs(:,1),EAs(:,2));
                                                                if in
                                                                    val = find(strcmpi("Eastern Asia",regionnames));
                                                                else
                                                                    in = inpolygon(lat,lon,Aus(:,1),Aus(:,2));
                                                                    if in
                                                                        val = find(strcmpi("Australia/New Zealand",regionnames));
                                                                    else
                                                                        error(strcat("Figure out where [",num2str(lat),",",num2str(lon),"] is"))
                                                                    end
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            regionnames(val)
        end
        cond = 1;
    end
end
outputvar(latitude,longitude) = val;