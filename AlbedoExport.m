% Export Albedo Atlas in Matlab/geotiff arrays for further analysis

% This script is used to export the albedo atlas from Gao et al., 2014, into a more readily readable
%   format for analysis on other computers. It is currently in an ENVI (bip) file format and read
%   is slow. Also, instead of saving files per months, they will be saved by biome types.

clearvars

% Get general inputs
computer = "GEO-005199";                % where the routine is run (for directory paths):
switch computer                         % directory paths for project data (depending on computer)
    case "GEO-007264"
        rootdir = 'C:\Work\GlobalAlbedo\';
        resultdir = rootdir;
    otherwise
        rootdir = 'D:\NCS-GlobalAlbedo\FilledAlbedo\';
        resultdir = 'E:\NCS-GlobalAlbedo\';
end

monthsnames = ["January","February","March","April","May","June","July","August","September",...
    "October","November","December"];
IGBPnn = ["Water","Evergreen_Needleleaf","Evergreen_Broadleaf","Deciduous_Needleleaf",...
    "Deciduous_Broadleaf","Mixed_Forests","Closed_Shrublands","Open_Shrublands","Woody_Savannas",...
    "savannas","Grasslands","Wetlands","Croplands","Urban","Crop_NaturalVeg_Mosaic","Snow_Ice",...
    "Barren"];




glodatafname = strcat(rootdir,"AlbedoGeneralData.mat");
load(glodatafname,'regoutputfiles','nblocks','nlon','nlat','IGBPBiomes','biomes_igbp',...
    'nbiomes','nmonths','albsn','nalb','presland','albedo_types','nblx','blocksize')
% save(,,'boxpath','resultslocalfolder',...
%     'nblx','nbly','nkernels','Ca_PgC',...
%     'nfrst','monthnames','albid','fillvalue_albedo','scale_albedo',...
%     'kernels','kernelscale','Aglobe','latlonscale',...
%     'pathways','pwnames','npw','pathwayslandcover','nlcc','pwid');

% Check that my new IGBP names are in the same sequence as the original ones
for ff = 1 : nbiomes
    nn = strsplit(IGBPnn(ff),"_");
    for k = 1 : length(nn)
        if nn(k) == "NaturalVeg", nn(k) = "natural vegetation"; end
        i = contains(IGBPBiomes(ff),nn(k),'IgnoreCase',true);
        if i == false
            error("My new IGBP names do not correspond to the ones in the original files")
        end
    end
end


albedofolder = strcat(resultdir,"AlbedoMaps\");
if exist(albedofolder,'dir') ~= 7
    mkdir(albedofolder)
end
albedomaps = NaN(nlat,nlon,nmonths);
R = georasterref('RasterSize',[nlat,nlon],'RasterInterpretation','cells','ColumnsStartFrom',...
    'north','LatitudeLimits',[-90,90],'LongitudeLimits',[-180,180]);
cmptag.Compression = 'LZW';             % LZW compression is better than the default ...    


% Read and write data (I'll loop through albedo types first for memory management)
for aa = 1 : nalb
    albedoname = albsn(aa);
    albedofilenames = strings(nbiomes,1);
    for ff = 1 : nbiomes
        fname = strcat(biomes_igbp(ff),"_",albedoname);
        albedofilenames(ff) = fname;
        eval(strcat(fname," = albedomaps;"));
    end
        
    for rr = 1 : nblocks
        if presland(rr) == true
            subdatafname = strcat(regoutputfiles,"AlbedoSubData_",num2str(rr),".mat");
            load(subdatafname,'regi','regj','blocklandmask',albedoname)
            eval(strcat("albedo = ",albedoname,";"))
            
            for ff = 1 : nbiomes
                bioname = biomes_igbp(ff);
                switch bioname
                    case "WAT"
                        if rr == 1
                            oceanvalues = NaN(nmonths,1);
                            for mm = 1 : nmonths
                                alb = albedo(:,:,mm,ff);
                                val = unique(alb(blocklandmask==0));
                                if numel(val) == 1
                                    oceanvalues(mm) = val;
                                else
                                    val2 = val(isfinite(val));
                                    if numel(val2) == 1
                                        oceanvalues(mm) = val2;
                                    else
                                        map = alb(blocklandmask==0);
                                        x = histcounts(map,[val2;1]);
                                        val3 = val2(x==max(x));
                                        if numel(val3) == 1
                                            oceanvalues(mm) = val3;
                                        else
                                        error(strcat("I have more than one value for oceans ",...
                                            "in ",albedo_types(aa)," albedos in ",monthsnames(mm)))
                                        end
                                    end
                                end
                            end
                        end                        
                    otherwise
                        for mm = 1 : nmonths
                            alb = albedo(:,:,mm,ff);
                            alb(blocklandmask==0) = NaN;
                            albedo(:,:,mm,ff) = alb; %#ok<SAGROW>
                        end
                end
                eval(strcat(bioname,"_",albedoname,"(regi,regj,:) = albedo(:,:,:,ff);"));
            end
        else
            ff = 1;
            y = ceil(rr/nblx);
            x = rr - (y-1)*nblx;
            regj = (x-1)*blocksize + 1 : x * blocksize;
            regi = (y-1)*blocksize + 1 : y * blocksize;
            albedo = repmat(permute(oceanvalues,[3,2,1]),[blocksize,blocksize,1]);
            eval(strcat(biomes_igbp(ff),"_",albedoname,"(regi,regj,:) = albedo;"));
            strcat("Should have done the albedo calculation for oceans in non-land block ",...
                num2str(rr)," (",biomes_igbp(ff),")")
        end
        strcat("Done for ",albedo_types(aa)," subregion ",num2str(rr))
    end
    
    % Save in Matfile
    matfile = strcat(resultdir,"AlbedoAtlas.",albedoname,".mat");
    save(matfile,'nlon','nlat','IGBPBiomes','biomes_igbp','nbiomes','nmonths',...
        albedofilenames{:},'-v7.3')
    
    % Save as geotiffs
    for ff = 1 : nbiomes
        abf = albedofilenames(ff);
        foldername = strcat(albedofolder,IGBPnn(ff),"\");
        if exist(foldername,'dir') ~= 7
            mkdir(foldername)
        end
        
        for mm = 1 : nmonths
            mne = monthsnames(mm);
            fname = strcat(foldername,abf,"_",mne,".tif");
            eval(strcat("data = ",abf,"(:,:,mm);"))
            geotiffwrite(fname,data,R,'TiffTags',cmptag);
        end
    end
end


            
            

