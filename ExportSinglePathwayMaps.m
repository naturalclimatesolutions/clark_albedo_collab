% Export GeoTIFF maps of all single land-cover change pathway albedo-change induced CO2e TOA RF


[~,lats] = meshgrid(lon05,lat05);
latA = -60;
Antarctica = lats <= latA;


for pp = 1 : npw
    pname = pwnames(pp);
    ii = find(pwid==pp);
    fold = strcat(resultslocalfolder,pname,"\ByIndividualKernel\");
    if ~exist(fold,'dir')
        mkdir(fold)
    end
    fold2 = strcat(resultslocalfolder,pname,"\ByStats\");
    if ~exist(fold2,'dir')
        mkdir(fold2)
    end
    for ll = 1 : length(ii)
        z = ii(ll);
        for kk = 1 : nkernels
            shortname = strcat(pathwayslandcover(z,1),"2",pathwayslandcover(z,2),"_",...
                kernels(kk),".tif");
            eval(strcat("data = ",pname,"(:,:,kk,ll);"))
            data(Antarctica) = NaN; %#ok<SAGROW>
            data(landmask==0) = NaN; %#ok<SAGROW> % this should not be necessary, but somehow I still got zero values
            fname = strcat(fold,shortname);
            geotiffwrite(fname,data,R05,'TiffTags',cmptag);
            strcat("Done with exporting ",shortname)
        end
        for ss = 1 : nstats
            eval(strcat("data =",pname,"Kstats(:,:,ss,ll);"))
            des = desiredstatistics(ss);
            des = replace(des," ","_");
            data(Antarctica) = nan;
            data(landmask==0) = nan;
            shortname = strcat(pathwayslandcover(z,1),"2",pathwayslandcover(z,2),"_",des,".tif");
            fname = strcat(fold2,shortname);
            geotiffwrite(fname,data,R05,'TiffTags',cmptag);
            strcat("Done with exporting ",shortname)
        end
    end
    clear pname ii ll z data fname fold fold2 shortname
end