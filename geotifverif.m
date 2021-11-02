% Compare exported geotiff data with original data


clearvars

computer = "GEO-005199";                % where the routine is run (for directory paths):

switch computer
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
end

glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'IGBPBiomes','biomes_igbp','latlonscale','nblocks','nlat','nlon','nlcc',...
    'nmonths','nkernels','npw','presland','pathwayslandcover','pwnames','pwid','kernels',...
    'regoutputfiles','resultslocalfolder','landmask','blocksize','kernelscale')

datafname = strcat(regoutputfiles,"AlbedoFinalData6ker.mat");

lat05 = 90 - latlonscale/2 : -latlonscale : -90 + latlonscale/2;
lon05 = -180 + latlonscale/2 : latlonscale : 180 - latlonscale/2;
[lons,lats] = meshgrid(lon05,lat05);

greyval = [0.7,0.7,0.7];
lightgreyval = [.9,.9,.9];
hotmap = colormap(hot);
jetmap = colormap(jet);
lastblue = find(jetmap(:,1)==0,1,'last');
firstyellow = find(jetmap(:,1)==1,1,'first');

valcat = [-1200,-800:200:-200,-100,-50,-20,20,50,100,200:200:800,1200];
nbcc = floor((length(valcat)-1)/2);
blueint = floor(lastblue/nbcc);
bluescale = jetmap(1:blueint:lastblue,:);
nlb = numel(bluescale)/3;
if nlb > nbcc, bluescale = bluescale([1:nlb-2,nlb],:); end
yellowint = floor((length(jetmap)-firstyellow)/nbcc);
yellowscale = jetmap(firstyellow:yellowint:length(jetmap),:);
nly = numel(yellowscale)/3;
if nly > nbcc, yellowscale = yellowscale([1,3:nly],:); end
catblind = cat(1,bluescale,greyval,yellowscale);



for ll = 1 : nlcc
    name = strcat(pathwayslandcover(ll,1),"2",pathwayslandcover(ll,2));
    lc1 = biomes_igbp == pathwayslandcover(ll,1);
    lc2 = biomes_igbp == pathwayslandcover(ll,2);

    folder = pwnames(pwid(ll));
    
    for kk = 1 : nkernels
        tifdata = readgeoraster(strcat(resultslocalfolder,folder,"\",name,"_",kernels(kk),".tif"));
        intline = find(find(pwid == pwid(ll)) == ll);
        
        if exist(folder,'var') == 0
            load(datafname,folder)
        end
        eval(strcat("matdata = ",folder,"(:,:,kk,intline);"));
        
        difmap = matdata - tifdata;
        dum = difmap;
        dum(isnan(dum)) = 0;
        
        if isempty(find(dum,1)), continue; else, warning("I have differences in my data!"); end
        
        catdata = discretize(difmap,valcat);
        
        val1 = round(min(difmap,[],'all'));
        val2 = round(max(difmap,[],'all'));
        
        
        figure(10*kk+ll); clf
        h.Units = 'pixels';
        h.Position = [1,31,2560,1333];
        axesm('MapProjection','robinson','Frame','on','MeridianLabel','on','ParallelLabel','on');
        load coast
        pcolorm(lats,lons,catdata)
        caxis([1 length(valcat)])
        colormap(catblind)
        c = colorbar;
        c.Ticks = 1:length(valcat);
        c.TickLabels = num2cell(valcat);
        c.FontSize = 20;
        c.Label.String = 'CO_2e [t CO_2 ha^-^1]';
        c.Ticks = 1:length(valcat);
        c.TickLabels = num2cell(valcat);
        plotm(lat,long,'Color','Black')
        ax = gca;
        ax.Title.FontSize = 24;
        ax.Title.String = {strcat(IGBPBiomes(lc1)," to ",IGBPBiomes(lc2));...
            strcat("\rm CO_2e differences between original and exported data (",...
            kernels(kk),") - range is [",num2str(round(val1)),",",num2str(round(val2)),"]")};
        gridm
        fname = strcat(resultslocalfolder,"Figures\",pathwayslandcover(ll,1),"2",...
            pathwayslandcover(ll,2),".",kernels(kk),".differences");
        print(strcat(fname,".jpg"),"-djpeg")
        
        pause
    end
    
end