function CO2eq_TOA = ptradforcing(initiallandcover,finallandcover,...
    nosnowblacksky,nosnowwhitesky,snowcovblacksky,snowcovwhitesky,...
    snowcover,diffuseradiation,kernelvalue,pixelarea,kernelstat)

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
    
ta = zeros(nmonths,nker);
Aglobe  = 5.1007e14;                    % global surface area [m2]
Ca_PgC  = 400*2.13;                     % current base CO2 in atmos [Pg C]


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
    
    if snow == 0
        Wda = squeeze(daNSBS .* (1 - diffuse) + daNSWS .* diffuse);
    elseif snow == 1
        Wda = squeeze(daSCBS .* (1 - diffuse)  + daSCWS .* diffuse);
    else
        Wda = squeeze(...
            daNSBS .* (1 - snow) .* (1 - diffuse) + ...
            daNSWS .* (1 - snow) .* diffuse + ...
            daSCBS .* snow .* (1 - diffuse) + ...
            daSCWS .* snow .* diffuse);
    end
    
    ta(m,:) =  Wda .* kernelvalue(m,:);        
    
end

arf = squeeze(mean(ta,1,'omitnan'));

if numel(arf) > 1
    sarf = zeros(length(kernelstat),1);
    for ss = 1 : length(kernelstat)
        eval("kstat = kernelstat(ss);")
        eval(strcat("sarf(ss) = ",kstat,"(arf);"))
    end
else
    sarf = arf;
end

% b. Global annual mean RF
glob = sarf .* pixelarea ./ Aglobe;

% c. CO2e (Pg C) (Change sign to make cooling a positive value)
PgC_TOA = Ca_PgC .* exp(-1 .* glob ./ 5.35);

% d. in tonnes CO2 ha-2
CO2eq_TOA = (PgC_TOA - Ca_PgC) .* 1e13 ./ pixelarea .* 44/12;
% CO2_2 = 1/0.89 .* arf .*10; checked that it is very similar to the other calculation...


