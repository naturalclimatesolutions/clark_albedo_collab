% Export GeoTIFF maps of all single land-cover change pathway albedo-change induced CO2e TOA RF

% Prepare arrays
[~,lats] = meshgrid(lon05,lat05);
latA = -60;
Antarctica = lats <= latA;
emptymap = ones(nlat05,nlon05,'int16') .* i16noval;
statsnames = replace(singlepathstatistics," ","");

for ll = 1 : nlcc
    for kk = 1 : nkernels
        kmapname = strcat(lower(pathwayslandcover(ll,1)),lower(pathwayslandcover(ll,2)),kernels(kk));
        eval(strcat(kmapname," = emptymap;"))
    end
    for ss = 1 : nstats
        smapname = strcat(lower(pathwayslandcover(ll,1)),lower(pathwayslandcover(ll,2)),statsnames(ss));
        eval(strcat(smapname," = emptymap;"))
    end
end

% Create world maps
for rr = 1 : nblocks
    if presland(rr) == false, continue; end
    subdatafname = strcat(regoutputfiles,"Albedo05_",num2str(rr),".mat");
    load(subdatafname,"regco2e","regco2estat","regi","regj")

    for ll = 1 : nlcc
        for kk = 1 : nkernels
            kmapname = strcat(lower(pathwayslandcover(ll,1)),lower(pathwayslandcover(ll,2)),kernels(kk));
            eval(strcat(kmapname,"(regi,regj) = int16(regco2e(:,:,kk,ll));"))
        end
        for ss = 1 : nstats
            smapname = strcat(lower(pathwayslandcover(ll,1)),lower(pathwayslandcover(ll,2)),statsnames(ss));
            eval(strcat(smapname,"(regi,regj) = int16(regco2estat(:,:,ss,ll));"))
        end
        clear kmapname smapname kk ss
    end

    if floor(rr/20) == rr/20
        strcat("done with building single pathway maps in subregion ",num2str(rr))
    end
end

% Export GeoTIFFs
for ll = 1 : nlcc
    pp = str2double(pathwayslandcover(ll,3));
    pname = pwnames(pp);
    fold = strcat(resultsGlasgowCOP,pname,"\ByIndividualKernel\");
    if ~exist(fold,'dir')
        mkdir(fold)
    end
    fold2 = strcat(resultsGlasgowCOP,pname,"\ByStats\");
    if ~exist(fold2,'dir')
        mkdir(fold2)
    end
    arootname = strcat(lower(pathwayslandcover(ll,1)),lower(pathwayslandcover(ll,2)));
    grootname = strcat(pathwayslandcover(ll,1),"2",pathwayslandcover(ll,2));

    for kk = 1 : nkernels
        kmapname = strcat(arootname,kernels(kk));
        kgeofiffname = strcat(fold,grootname,"_",kernels(kk),".tif");
        eval(strcat("data =",kmapname,";"))
        data(Antarctica) = i16noval; %#ok<SAGROW> 
        geotiffwrite(kgeofiffname,data,R05,'TiffTags',cmptag);
        clear kmapname kgeofiffname data
    end
    for ss = 1 : nstats
        smapname = strcat(arootname,statsnames(ss));
        sgeofiffname = strcat(fold2,grootname,"_",statsnames(ss),".tif");
        eval(strcat("data =",smapname,";"))
        data(Antarctica) = i16noval;  
        geotiffwrite(sgeofiffname,data,R05,'TiffTags',cmptag);
        clear smapname sgeofiffname data
    end
    strcat("done with exporting geotiff maps of ",arootname)
    clear pp pname fold fold2 arootname grootname kk ss
end




% Tifftags can be
% 'Artist', 'Compression', 'Copyright', 'DateTime', 'DocumentName', 'DotRange', 'ExtraSamples', 'FillOrder', 'Group3Options',
% 'Group4Options', 'HalfToneHints', 'HostComputer', 'ICCProfile', 'ImageDepth', 'ImageDescription', 'InkNames', 'InkSet',
% 'JPEGColorMode', 'JPEGQuality', 'Make', 'MaxSampleValue', 'MinSampleValue', 'Model', 'NumberOfInks', 'Orientation', 'PageName',
% 'PageNumber', 'PhotometricInterpretation', 'Photoshop', 'PlanarConfiguration', 'PrimaryChromaticities', 'RPCCoefficientTag',
% 'ReferenceBlackWhite', 'ResolutionUnit', 'RichTIFFIPTC', 'RowsPerStrip', 'SGILogDataFmt', 'SMaxSampleValue', 'SMinSampleValue',
% 'SToNits', 'Software', 'TargetPrinter', 'Thresholding', 'TileLength', 'TileWidth', 'TransferFunction', 'WhitePoint', 'XMP',
% 'XPosition', 'XResolution', 'YCbCrCoefficients', 'YCbCrPositioning', 'YCbCrSubSampling', 'YPosition', 'YResolution', 'ZipQuality'
