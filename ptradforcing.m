function [CO2eq_TOA,GWPRF] = ptradforcing(initiallandcover,finallandcover,...
    nosnowblacksky,nosnowwhitesky,snowcovblacksky,snowcovwhitesky,...
    snowcover,diffuseradiation,kernelvalue,pixelarea,kernelstat,gwp100,capgc)

% FUNCTION ptradforcing: rafiative forcing at a particular pixel
%   This function calculate the annual mean albedo-induced radiative forcing CO2 equivalent at
%   top-of-atmosphere for a given pixel. All albedo, snow cover, diffuse radiation fractions and
%   kernel values need to be provided as monthly arrays. If more than one kernel is provided, then
%   one or several statistic should also be defined, if not, median is the default.


[nmonths,nker] = size(kernelvalue);
if nmonths ~= 12
    if nker == 12
        kernelvalue = permute(kernelvalue,[2,1]);
        [nmonths,nker] = size(kernelvalue);
    else
        error("I have an unexpected array size for kernelvalue!")
    end
end

if nargin < 11
    kernelstat = "median";
end
if nargin < 12
    gwp100 = 1.9;
end
if nargin < 13
    capgc = 400*2.13;                   % current base CO2 in atmos [Pg C]
end
    
ta = zeros(nmonths,nker);
Aglobe  = 5.1007e14;                    % global surface area [m2]


for m = 1 : nmonths
    
    if finallandcover == 22             % mean between grass and crop (asked as layer by TNC)
        finnsbs = mean([nosnowblacksky(m,11),nosnowblacksky(m,13)]);
        finnsws = mean([nosnowwhitesky(m,11),nosnowwhitesky(m,13)]);
        finscbs = mean([snowcovblacksky(m,11),snowcovblacksky(m,13)]);
        finscws = mean([snowcovwhitesky(m,11),snowcovwhitesky(m,13)]);
    else
        finnsbs = nosnowblacksky(m,finallandcover);
        finnsws = nosnowwhitesky(m,finallandcover);
        finscbs = snowcovblacksky(m,finallandcover);
        finscws = snowcovwhitesky(m,finallandcover);
    end
    
    daNSBS = finnsbs- nosnowblacksky(m,initiallandcover);
    daNSWS = finnsws - nosnowwhitesky(m,initiallandcover);
    daSCBS = finscbs - snowcovblacksky(m,initiallandcover);
    daSCWS = finscws - snowcovwhitesky(m,initiallandcover);
    
    snow = snowcover(m);
    diffuse = diffuseradiation(m);

    Wda = squeeze(...
        daNSBS .* (1 - snow) .* (1 - diffuse) + daNSWS .* (1 - snow) .* diffuse + ...
        daSCBS .* snow .* (1 - diffuse) + daSCWS .* snow .* diffuse);
        
    ta(m,:) =  Wda .* kernelvalue(m,:);
    
end

arf = squeeze(mean(ta,1,'omitnan'));
glob = arf .* pixelarea ./ Aglobe;          % Global annual mean RF
PgC_TOA = capgc .* exp(-1 .* glob ./ 5.35);% C02e (Pg C) (Change sign to make cooling positive)
co2 = (PgC_TOA - capgc) .* 1e13 ./ pixelarea .* 44/12; % in MgCo2/ha
% CO2_2 = 1/0.89 .* arf .*10; checked that it is very similar to the other calculation...

% Calculate kernel statistics
if nker > 1 && sum(isnan(co2)) < nker
    CO2eq_TOA = zeros(length(kernelstat),1);
    vsig = sign(co2);
    if sum(vsig) == numel(co2)
        minval = min(co2); maxval = max(co2);
    else
        maxval = min(co2);
        mi = min(abs(co2));
        ii = abs(co2) == mi;
        tmpminval = co2(ii);
        if numel(tmpminval) == 1
            minval = tmpminval;
        elseif sum(diff(tmpminval)) == 0
            minval = tmpminval(1);
        else
            minval = mi;
        end
    end
    for ss = 1 : length(kernelstat)
        kstat = kernelstat(ss);
        switch kstat
            case "min"
                CO2eq_TOA(ss) = minval;
            case "max"
                CO2eq_TOA(ss) = maxval;
            otherwise
                eval(strcat("CO2eq_TOA(ss) = ",kstat,"(co2);"))
        end
    end
elseif sum(isnan(co2)) == nker
    CO2eq_TOA = nan;
else
    CO2eq_TOA = co2;
end

GWPRF = CO2eq_TOA .* gwp100;


