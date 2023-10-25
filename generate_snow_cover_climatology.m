% Generate Global Snow Cover Climatology from MODIS data
%
% Generage climatology by averaging all available yearly data from
%  MODIS/Terra Snow Cover Monthly L3 Global 0.05Deg CMG, Version 6
%   This data set reports monthly average snow cover in 0.05 degree (approx. 5 km) resolution 
%   Climate Modeling Grid (CMG) cells. Monthly averages are computed from daily snow cover 
%   observations in the MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG (MOD10C1) data set.
%   
%  Source: https://nsidc.org/data/modis/data_summaries
%  Data Set ID: MOD10CM (version 6.1)
%
% Created by Prof. Chris Williams, Clark University


% Define input variables and parameters
% *************************************
daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31];
firstdayinmonth = cumsum([1,daysinmonth]);
lydaysinmonth = daysinmonth;
lydaysinmonth(2)= 29;
lyfirstdayinmonth = cumsum([1,lydaysinmonth]);
daymonthindex = cat(1,firstdayinmonth,lyfirstdayinmonth);
latlonscale = 0.05;
nlat = 180 / latlonscale;
nlon = 360 / latlonscale;
blankmap = zeros(nlat,nlon);
clear daysinmonth firstdayinmonth lydaysinmonth lyfirstdayinmonth

if ~exist(climsnowcover,'dir'), mkdir(climsnowcover); end
filelist = deblank(string(ls(snowcoverdata)));
snowfilelist = filelist(endsWith(filelist,"hdf"));

yeardata = 2000:2022;   % add years when more available ...
availablefiles = strings(length(yeardata),length(monthnames));


% 2. Loop over months
% -------------------
for m = 1 : length(monthnames)
    
    data = blankmap;
    validdata = blankmap;
    dayindex = unique(daymonthindex(:,m));

    
    % 3. Loop through all files and see if they belong to the desired month
    % ---------------------------------------------------------------------
    for n = 1 : length(snowfilelist)
        
        fname = snowfilelist(n);
        date = char(extractBetween(fname,"MOD10CM.A",".061"));
        fyear = str2double(date(1:4));
        fday = str2double(date(5:7));
        if ismember(fday,dayindex) == 0, continue;end
        yi = fyear == yeardata;
        availablefiles(yi,m) = fname;
        
        
        % 4. Get data from file (if it belongs in desired month)
        % ---------------------
        S = hdfinfo(strcat(snowcoverdata,fname),'eos');
        varname = S.Grid.DataFields(1).Name;
        A = hdfread(S.Filename,varname);
        
        
        % 5. Check validity and add to previous data
        % ------------------------------------------
        A = double(A);
        A(A>100) = NaN;
        validpixels = isfinite(A);
        
        data = sum(cat(3,data,A),3,'omitnan');
        validdata = sum(cat(3,validdata,validpixels),3);
        
    end
    
    
    % 6. Calculate mean and save
    % --------------------------
    snowcover = data ./ validdata;
    
    outfname = strcat(climsnowcover,"SnowCover_2000-2021_",monthnames(m),".mat");
    save(outfname,'snowcover')
    
end

vartype = repmat("string",[1,length(monthnames)+1]);
T = table('Size',[length(yeardata),length(monthnames)+1],'VariableTypes',vartype);
T(:,1) = cellstr(string(yeardata)');
for ii = 2 : length(monthnames)+1
    T(:,ii) = cellstr(availablefiles(:,ii-1));
end
T.Properties.RowNames = cellstr(string(yeardata)); % This is redundant, but rowname is not saved in excel
T.Properties.VariableNames = cellstr(["year",monthnames]);
writetable(T,strcat(climsnowcover,"SnowFileNames.xlsx"));

