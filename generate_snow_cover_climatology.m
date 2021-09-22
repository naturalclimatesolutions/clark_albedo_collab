% Generate Global Snow Cover Climatology from MODIS data
%
% Generage climatology by averaging all available yearly data from
%  MODIS/Terra Snow Cover Monthly L3 Global 0.05Deg CMG, Version 6
%   This data set reports monthly average snow cover in 0.05 degree (approx. 5 km) resolution 
%   Climate Modeling Grid (CMG) cells. Monthly averages are computed from daily snow cover 
%   observations in the MODIS/Terra Snow Cover Daily L3 Global 0.05Deg CMG (MOD10C1) data set.
%   
%  Source: https://nsidc.org/data/MOD10CM/versions/6
%  Data Set ID: MOD10CM
%
% Created by Prof. Chris Williams, Clark University

% Last modified: Natalia Hasler 9/22/2021
%                (See Modification history at end of file)



clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1. Define input variables and parameters
% ----------------------------------------

monthnames = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"];
daysinmonth = [31,28,31,30,31,30,31,31,30,31,30,31];
firstdayinmonth = cumsum([1,daysinmonth]);
lydaysinmonth = daysinmonth;
lydaysinmonth(2)= 29;
lyfirstdayinmonth = cumsum([1,lydaysinmonth]);
latlonscale = 0.05;
nlat = 180 / latlonscale;
nlon = 360 / latlonscale;
blankmap = zeros(nlat,nlon);

MODISSnowFilesFolderPath = "G:TNC\GlobalAlbedo\SnowData\MODIS_SNOW_DATA\";
ResultFolderPath = "G:TNC\GlobalAlbedo\SnowData\";
filelist = deblank(string(ls(MODISSnowFilesFolderPath)));
snowfilelist = filelist(endsWith(filelist,"hdf"));


% 2. Loop over months
% -------------------
for m = 1 : length(monthnames)
    
    data = blankmap;
    validdata = blankmap;

    
    % 3. Loop through all files and see if they belong to the desired month
    % ---------------------------------------------------------------------
    for n = 1 : length(snowfilelist)
        
        fname = snowfilelist(n);
        a = replace(fname,"."," ");
        fdate = erase(extractBetween(a,"MOD10CM "," 006 "),"A");
        fyear = floor(str2double(fdate) / 10^3);
        fday = str2double(fdate) - fyear * 10^3;
        
        if isLeapYear(fyear)
            day1 = lyfirstdayinmonth(m); day2 = lyfirstdayinmonth(m+1);
        else
            day1 = firstdayinmonth(m); day2 = firstdayinmonth(m+1);
        end
        
        if fday < day1 || fday >= day2; continue; end
        
        
        % 4. Get data from file (if it belongs in desired month)
        % ---------------------
        S = hdfinfo(strcat(MODISSnowFilesFolderPath,fname),'eos');
        varname = S.Grid.DataFields.Name;
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
    
    outfname = strcat(ResultFolderPath,"SnowCover_2000-2021_",monthnames(m),".mat");
    save(outfname,'snowcover')
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODIFICATIONS

% 9/2/21 NH Corrected the calculation of valid data
