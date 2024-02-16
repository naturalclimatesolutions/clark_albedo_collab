function [area,co2,seq,missingareas,mismask] = binnedvalues(input,carbon,nci,mask,categories,...
    negco2,region,ecoregion,pixelarea,regionlist,ecoregionlist,...
    regionnodata,ecoregionnodata,datamask,areascalar,co2scalar,inputscalar)

cllist = ["uint8","int8","uint16","int16","uint32","int32","single","double","logical"];
mislist = [2^8-1,-2^7,2^16-1,-2^15,2^32-1,-2^31,nan,nan,false];

if nargin < 17
    inputscalar = 10^4;
end
if nargin < 16
    co2scalar = 10^9;
end
if nargin < 15
    areascalar = inputscalar * 10^6;
end
if nargin < 14
    datamask = true(size(input));
end
if nargin < 13
    ecoregionnodata = mislist(cllist == class(ecoregion));
end
if nargin < 12
    regionnodata = mislist(cllist == class(region));
end

if isa(negco2,'logical')
    negco2 = find(negco2);
end

area = zeros(numel(ecoregionlist),numel(categories)-1,numel(regionlist)+1);
co2 = zeros(numel(ecoregionlist),numel(categories)-1,numel(regionlist)+1);
seq = zeros(numel(ecoregionlist),numel(categories)-1,numel(regionlist)+1);
missingareas = zeros(numel(ecoregionlist)+1,numel(regionlist)+1);
mismask = false(size(input));

thisblockeco = setdiff(unique(ecoregion(datamask)),ecoregionnodata);
thisblockreg = unique(region(datamask));
if ismember(0,thisblockreg)
    region(region==0) = regionnodata;
    thisblockreg = unique(region(datamask));
end

ormask = mask;
mask(logical(1-datamask)) = false;
difmask = logical(ormask - mask);
difmaskarea = sum(pixelarea(difmask));

input(input<categories(1)) = categories(1);
input(input>categories(numel(categories))) = categories(numel(categories));
catinput = discretize(input,categories);

verifpix = zeros(numel(thisblockeco),numel(thisblockreg));

for ee = 1 : numel(thisblockeco)
    eco = thisblockeco(ee);
    ecopix = ecoregion == eco;
    ecoi = ecoregionlist==eco;
    for rr = 1 : numel(thisblockreg)
        reg = thisblockreg(rr);
        regpix = region == reg;
        if reg == regionnodata, regi = numel(regionlist)+1; else, regi = regionlist==reg; end

        selpix = ecopix + regpix + mask == 3;

        verifpix(ee,rr) = sum(selpix,'all');
        selectedareas = pixelarea(selpix);

        selectinput = catinput(selpix);

        no = histcounts(selectinput,1:numel(categories)); nok = find(no);
        if sum(no) ~= sum(selpix,'all')
            kk = isnan(selectinput);
            if sum(kk,'all') == sum(selpix,'all') - sum(no)
                missingareas(ecoi,regi) = sum(selectedareas(kk),'all') ./ areascalar;
                mm = isnan(input); zz = mm + selpix == 2;
                mismask(zz) = true;
            else
                error("figure this one out: not nan input")
            end
        end
        
        selectcarbon = carbon(selpix);
        selectnci = nci(selpix);

        for ll = 1 : numel(nok)
            clo = nok(ll);
            pts = selectinput == clo;
            if sum(pts) ~= no(clo), error("My count is wrong in albedo"); end
            area(ecoi,clo,regi) = sum(selectedareas(pts),'omitnan') ./ areascalar;
            binco2 = sum(selectnci(pts) .* selectedareas(pts) ./ inputscalar,'omitnan') ./co2scalar;
            if ismember(clo,negco2) && binco2 > 0 %#ok<BDSCI>
                error("My carbon total is not correct")
            end 
            co2(ecoi,clo,regi) = binco2;
            seq(ecoi,clo,regi) = sum(selectcarbon(pts) .* selectedareas(pts) ./ inputscalar,'omitnan') ./co2scalar;
        end
    end

    if sum(verifpix(ee,:)) == sum(ecopix(mask))
            area(ecoi,:,numel(regionlist)+1) = sum(area(ecoi,:,:),3);
            co2(ecoi,:,numel(regionlist)+1) = sum(co2(ecoi,:,:),3);
            seq(ecoi,:,numel(regionlist)+1) = sum(seq(ecoi,:,:),3);
    else
        error("figure this one out : verifpix in ecoregion")
    end
end

if sum(verifpix,'all') ~= sum(mask,'all')
    nbpix = sum(mask,'all') - sum(verifpix,'all');
    ndp = (ecoregion==ecoregionnodata) + mask == 2;
    if sum(ndp,'all') == nbpix
        tte = unique(ecoregion(ndp));
        ttr = unique(region(ndp));
        pcs = 0;
        for ee = 1 : numel(tte)
            eco = tte(ee);
            if eco == ecoregionnodata, ecoi = numel(ecoregionlist)+1;
            else, ecoi = ecoregionlist == eco; end
            for rr = 1 : numel(ttr)
                reg = ttr(rr);
                if reg == regionnodata, regi = numel(regionlist)+1;
                else, regi = regionlist(reg); end
                pixs = (ecoregion==eco) + (region==reg) + mask == 3;
                missingareas(ecoi,regi) = sum(pixelarea(pixs),'all') ./ areascalar;
                pcs = pcs + sum(pixs,'all');
                mismask(pixs) = true;
            end
        end
        if pcs ~= nbpix, error("recount the missing pixels in missing ecoregion/region"); end
    else
        error("Figure this one out: verifpix globally")
    end
end

missingareas(numel(ecoregionlist)+1,numel(regionlist)+1) = ...
    missingareas(numel(ecoregionlist)+1,numel(regionlist)+1) + difmaskarea;






