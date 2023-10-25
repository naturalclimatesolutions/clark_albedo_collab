function h = albedofigure(data,categories,latitudes,longitudes,figurename,positions,prctvals,...
    coastlat,coastlon,axislegend,mapcolors,climatearrows,lonflag)

% ALBEDOFIGURE
%   This function is used to create a Robinson projection world map of "data". It has the option of
%   drawing an area percentile bar next to the colorbar as well as arrows indicating the direction
%   of climate warming/cooling for easier interpretation. If the desired figure is part of a
%   three-figure panel, longitude labels are omitted.
%
% Input variables are
%   - data              % variable in [lat/lon] to be plotted
%   - categories        % variable categories (bins) for the colorbar
%   - latitudes         % latitude coordinates of each pixel (2D array, i.e. nlat,nlon)
%   - longitudes        % longitude coordinates of each pixel (2D array, i.e. nlat,nlon)
%   - figurename        % figure name (for saving)
%   - positions         % (optional) axes and figure positions
%   - prctvals          % (optional) pre-calculated areas statistics at [5,10,25,50,75,90,95]
%                           percentiles. If not passed, this will be calculated on the pixels array
%                           and would differ from the actual area percentiles. Code should be
%                           amended if the real area needs to be calculated here
%   - coastlat          % (optional) latitude coordinates of world coastlines
%   - coastlon          % (optional) longitude coordinates of world coastlines
%   - axislegend        % colorbar unit. Default is Mg CO2e ha-1
%   - mapcolors         % (optional) colorscheme for colorbar
%   - climatearrows     % (logical) flag to include (true) or not (false) the arrows next to the
%                           colorbar to indicate climate warming/cooling direction
%   - lonflag           % (logical) flag to include (true) or exclude (false) longitude lables.
%                           Default is true unless the figure name includes "b_" or "c_"
%
% Created by Natalia Hasler, Clark University -  September 2022


if nargin < 13
    if contains(figurename,'b_') || contains(figurename,'c_')
        lonflag = false;
    else
        lonflag = true;
    end
end
if nargin < 12
    climatearrows = true;
end   
if nargin < 11
    mapcolors = redtobluecolorbar(categories);
end
if nargin < 10
    axislegend = 'Mg CO_2e ha^-^1';
end
if nargin < 8
    load coastlines %#ok<LOAD> 
end
if nargin < 6
    figpos = [1 31 2560 1333];
    axpos = [0.0426 0.0190 0.8344 0.9060];
else
    figpos = positions.figure;
    axpos = positions.axes;
end

cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double","logical"];
mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan,false];

dataclass = class(data);
nodatavalue = mislist(strcmp(cllist,dataclass));
data = double(data);
[l,c] = size(data); onedimdata = reshape(data,[l*c,1]);

data(data==nodatavalue) = NaN;
catdata = discretize(data,categories);
catlabels = num2cell(categories);
catlabels(1) = num2cell(strcat("< ",num2str(categories(2))));
catlabels(numel(categories)) = num2cell(strcat("> ",num2str(categories(numel(categories)-1))));

prctscale = [5,10,25,50,75,90,95];
if nargin < 7 || numel(prctvals)==0
    prctvals = prctile(onedimdata,prctscale);   % with lots of nan, it works better in vector form
end

sclimits = prctile(onedimdata,[1,99]);
fakecats = categories; fakecats(1) = sclimits(1); fakecats(numel(fakecats)) = sclimits(2);
if numel(prctvals) > 1
    prctflag = true;
    prctbins = discretize(prctvals,categories);
    inbinpos = prctbins;
    for ii = 1 : numel(prctvals)
        ind = prctbins(ii);
        if ismember(prctvals(ii),[categories(1),categories(numel(categories))])
            val = sign(prctvals(ii)) * inf;
            if ind == 1, inbinpos(ii) = 0.1; else, inbinpos(ii) = 0.9; end
            prctvals(ii) = val;
        else
            val = prctvals(ii);
            inbinpos(ii) = (val - fakecats(ind)) / (fakecats(ind+1) - fakecats(ind));
        end
    end
%     d = [1;diff(prctbins)]; b = unique(prctbins(d==0));
%     vertal = repmat("middle",[numel(prctvals),1]);
%     for k = 1 : numel(b)
%         bb = b(k);
%         kk = find(prctbins==bb);
%         vertal(kk(1)) = "top";
%         vertal(kk(numel(kk))) = "bottom";
%     end
else
    prctflag = false;
end



clf
h = gcf;
h.Units = 'pixels';
h.Position = figpos;
if lonflag
    axesm('MapProjection','robinson','Frame','on','MLineLocation',20,'PLineLocation',20,...
        'Grid','on','MeridianLabel','on','MLabelLocation',60,'ParallelLabel','on',...
        'PLabelLocation',20,'LabelRotation','on','FontSize',18);
else
    axesm('MapProjection','robinson','Frame','on','MLineLocation',20,'PLineLocation',20,...
        'Grid','on','ParallelLabel','on','PLabelLocation',20,'FontSize',18);
end

pcolorm(latitudes,longitudes,catdata)
plotm(coastlat,coastlon,'Color','Black')
caxis([1 length(categories)])
colormap(mapcolors)
c = colorbar;
c.Ticks = 1:length(categories);
c.TickLabels = catlabels;
c.FontSize = 25;
c.Label.String = axislegend;
ax = gca;
ax.Position = axpos;
cpos = c.Position;
int = cpos(4) / (numel(categories)-1);
if prctflag
    x = cpos(1) - ([2,1.5]*cpos(3));
    y = cpos(2) + (prctbins-1+inbinpos)*int;
    annotation('line',[x(2) x(2)],[y(1) y(numel(y))])
    for ii = 1 : numel(prctscale)
        if prctscale(ii) == 50, lw = 1.5; else, lw = 0.5; end
        annotation('textarrow',[x(1) x(2)],[y(ii) y(ii)],'String',...
            num2str(round(prctvals(ii))),'HeadStyle','none','FontSize',20,'LineWidth',lw)
%             'VerticalAlignment',vertal(ii))        
    end
    annotation('line',[x(2) x(2)],[y(prctscale==25) y(prctscale==75)],'LineWidth',1.5)
end
if climatearrows
    yclwm = [.21 .06]; yclcl = [.78 .93];
    ha = ["center","center"];
    zp = discretize(0,categories);
    if zp <= 2 || numel(categories)-zp <= 2
        zinpos = (0 - categories(zp)) / (categories(zp+1) - categories(zp));
        zpos = cpos(2) + (zp-1+zinpos)*int;
        mid = cpos(2) + (cpos(4)/2);
        ha = ["right","left"];
        if zpos > mid
            yclcl(1) = zpos;
            ltrd = cpos(2) + (cpos(4)/3);
            yclwm = [ltrd ltrd-.15];
        else
            yclwm(1) = zpos;
            utrd = cpos(2) + 2*cpos(4)/3;
            yclcl = [utrd utrd+.15];
        end
    end
    annotation('textarrow',[.95 .95],yclwm,'String',"climate warming",'FontSize',25,...
        'TextRotation',90,'VerticalAlignment','top','HorizontalAlignment',ha(1))
    annotation('textarrow',[.95 .95],yclcl,'String',"climate cooling",'FontSize',25,...
        'TextRotation',90,'VerticalAlignment','top','HorizontalAlignment',ha(2))
end
print(figurename,"-djpeg")


% Formulas for tick positions
% for middle points
% c.Ticks = minvalue+(extent-1)/(2*extent):(extent-1)/extent:maxvalue;
% for intervals
% c.Ticks = minvalue:(extent-1)/extent:maxvalue;


% Statistics with % labels
% if prctflag
%     cpos = c.Position;
%     x = cpos(1) - ([3.5,3,2.5]*cpos(3));
%     int = cpos(4) / (numel(categories)-1);
%     y = cpos(2) + (prctbins-1+inbinpos)*int;
% 
%     annotation('line',[x(2) x(2)],[y(1) y(numel(y))])
%     for ii = 1 : numel(prctscale)
%         annotation('textarrow',[x(1) x(2)],[y(ii) y(ii)],'String',...
%             num2str(round(prctvals(ii))),'HeadStyle','none','FontSize',16)
%         annotation('textarrow',[x(3) x(2)],[y(ii) y(ii)],'String',...
%             strcat(num2str(prctscale(ii)),"%"),'HeadStyle','none','FontSize',16)
%         
%     end
% end
