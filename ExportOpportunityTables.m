% Export Opportunity (and total biome) tables

% Export Tables with all Opportunities - base/low/high separate
below50 = offsetcat(1:nocat)<aothreshold;
above75 = offsetcat(2:nocat+1)>75;
above100 = offsetcat(2:nocat+1)>100;
semigeographicorder = [11,6,8,12,5,4,7,3,2,1,14,9,10,13];
opporder = ["Griscom","Bastin","Walker","TotalBiome"];
no = numel(opporder);
[~,noi] = ismember(opporder,allreforestationopp);
noi(noi==0) = 4;
xcelfile = strcat(resfolder,"ReforestationOpportunityStatsbybiomes.xlsx");
ncont = numel(continents);
clm = ["Biomes","Area","No-albedo carbon","Total NCI (carbon+albedo)",...
    "NCI average rate","Mean Albedo Offset",...
    "<50%AO Areas","<50%AO %Areas","<50%AO NCI","<50%AO NCI average rate",...
    ">75%AO Areas",">75%AO %Areas",">75%AO NCI",">75%AO NCI average rate",...
    ">100%AO Areas",">100%AO %Areas",">100%AO NCI",">100%AO NCI average rate"];
clu = ["","[Mha]","[Pg CO2e]","[Pg CO2e]","[Mg/ha]","[%]",...
    "[Mha]","[%]","[Pg CO2e]","[Mg/ha]","[Mha]","[%]","[Pg CO2e]","[Mg/ha]",...
    "[Mha]","[%]","[Pg CO2e]","[Mg/ha]"];

TA = zeros(nbiomes,numel(opporder),numel(bdnames));
TM = TA; PA = TA; PM = TA; NA = TA; NM = TA; A75 = TA; M75 = TA; TC = TA;

for ss = 1 : numel(bdnames)
	sensname = bdnames(ss);
    worldvalues = zeros(nbiomes*(no+1),numel(clm)-1);
    atbl = eval(strcat("AreasbyBiome_AO",sensname,"(semigeographicorder,:,noi,ncont+1);"));
    ctbl = eval(strcat("TotalCO2byBiome_AO",sensname,"(semigeographicorder,:,noi,ncont+1);"));
    stbl = eval(strcat("JustCarbonbyBiome_AO",sensname,"(semigeographicorder,:,noi,ncont+1);"));
    TA(:,:,ss) = squeeze(sum(atbl,2));
    TM(:,:,ss) = squeeze(sum(ctbl,2));
    PA(:,:,ss) = squeeze(sum(atbl(:,below50,:),2));
    PM(:,:,ss) = squeeze(sum(ctbl(:,below50,:),2));
    A75(:,:,ss) = squeeze(sum(atbl(:,above75,:),2));
    M75(:,:,ss) = squeeze(sum(ctbl(:,above75,:),2));
    NA(:,:,ss) = squeeze(sum(atbl(:,above100,:),2));
    NM(:,:,ss) = squeeze(sum(ctbl(:,above100,:),2));
    TC(:,:,ss) = squeeze(sum(stbl,2));
    
    for bb = 1 : nbiomes
        tar = TA(bb,:,ss);
        mit = TM(bb,:,ss);
        tco = TC(bb,:,ss);
        posar = PA(bb,:,ss);
        posmit = PM(bb,:,ss);
        a75 = A75(bb,:,ss);
        m75 = M75(bb,:,ss);
        negar = NA(bb,:,ss);
        negmit = NM(bb,:,ss);
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,1) = tar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,2) = tco;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,3) = mit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,4) = mit*co2scalar ./(tar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,5) = (tco-mit) ./ tco .* 100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,6) = posar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,7) = posar./tar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,8) = posmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,9) = posmit*co2scalar ./(posar*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,10) = a75;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,11) = a75./tar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,12) = m75;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,13) = m75*co2scalar ./(a75*areascalar) * perareascalar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,14) = negar;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,15) = negar./tar.*100;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,16) = negmit;
        worldvalues((bb-1)*(no+1)+2:(bb-1)*(no+1)+no+1,17) = negmit*co2scalar ./(negar*areascalar) * perareascalar;
    end
    fc = reshape(cat(1,biomenames(semigeographicorder)',repmat(opporder',[1,nbiomes])),[nbiomes*(no+1),1]);
    
    T = cell(nbiomes*(no+1)+2,numel(clm));
    T(1,:) = cellstr(clm);
    T(2,:) = cellstr(clu);
    T(3:nbiomes*(no+1)+2,1) = cellstr(fc);
    T(3:nbiomes*(no+1)+2,2:numel(clm)) = num2cell(worldvalues);
    
    TT = cell2table(T);
    writetable(TT,xcelfile,'Sheet',sensname)    
end


% Combined Table (rounding, no decimals)
% **************
rclm = [1:6,8:10,12,16];
nclm = clm(rclm);
nclu = clu(rclm);

TMR = TM.*co2scalar ./(TA*areascalar) * perareascalar;
TMO = (TC - TM) ./ TC .* 100;
RPA = PA ./ TA .* 100;
PMR = PM.*co2scalar ./(PA*areascalar) * perareascalar;
RA75 = A75 ./ TA .* 100;
RNA = NA ./ TA .* 100;
M75R = M75.*co2scalar ./(A75*areascalar) * perareascalar;
NMR = NM.*co2scalar ./(NA*areascalar) * perareascalar;

obioname = biomenames(semigeographicorder);

T = cell((nbiomes+1)*(no+1)+2,numel(nclm));
T(1,:) = cellstr(nclm);
T(2,:) = cellstr(nclu);
for bb = 1 : nbiomes
    T((bb-1)*(no+1)+3,1) = cellstr(obioname(bb));
    for oo = 1 : no
        T((bb-1)*(no+1)+3+oo,1) = cellstr(opporder(oo));
        T((bb-1)*(no+1)+3+oo,2) = num2cell(round(TA(bb,oo,2)));
        T((bb-1)*(no+1)+3+oo,3) = num2cell(round(TC(bb,oo,2)));
        T((bb-1)*(no+1)+3+oo,4) = cellstr(strcat(num2str(round(TM(bb,oo,2)),'%i')," [",...
            num2str(round(TM(bb,oo,3)),'%i')," - ",num2str(round(TM(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,5) = cellstr(strcat(num2str(round(TMR(bb,oo,2)),'%i')," [",...
            num2str(round(TMR(bb,oo,3)),'%i')," - ",num2str(round(TMR(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,6) = cellstr(strcat(num2str(round(TMO(bb,oo,2)),'%i')," [",...
            num2str(round(TMO(bb,oo,3)),'%i')," - ",num2str(round(TMO(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,7) = cellstr(strcat(num2str(round(RPA(bb,oo,2)),'%i')," [",...
            num2str(round(RPA(bb,oo,3)),'%i')," - ",num2str(round(RPA(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,8) = cellstr(strcat(num2str(round(PM(bb,oo,2)),'%i')," [",...
            num2str(round(PM(bb,oo,3)),'%i')," - ",num2str(round(PM(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,9) = cellstr(strcat(num2str(round(PMR(bb,oo,2)),'%i')," [",...
            num2str(round(PMR(bb,oo,3)),'%i')," - ",num2str(round(PMR(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,10) = cellstr(strcat(num2str(round(RA75(bb,oo,2)),'%i')," [",...
            num2str(round(RA75(bb,oo,3)),'%i')," - ",num2str(round(RA75(bb,oo,1)),'%i'),"]"));
        T((bb-1)*(no+1)+3+oo,11) = cellstr(strcat(num2str(round(RNA(bb,oo,2)),'%i')," [",...
            num2str(round(RNA(bb,oo,3)),'%i')," - ",num2str(round(RNA(bb,oo,1)),'%i'),"]"));
    end
end
T(nbiomes*(no+1)+3,1) = cellstr("Globally");
GTA = squeeze(sum(TA));
GTC = squeeze(sum(TC));
GTM = squeeze(sum(TM));
GPA = squeeze(sum(PA));
GPM = squeeze(sum(PM));
GA75 = squeeze(sum(A75));
GNA = squeeze(sum(NA));
GM75 = squeeze(sum(M75));
GNM = squeeze(sum(NM));
GTMR = GTM.*co2scalar ./(GTA*areascalar) * perareascalar;
GTMO = (GTC - GTM) ./ GTC .* 100;
GRPA = GPA ./ GTA .* 100;
GPMR = GPM.*co2scalar ./(GPA*areascalar) * perareascalar;
GRA75 = GA75 ./ GTA .* 100;
GRNA = GNA ./ GTA .* 100;
GM75R = GM75.*co2scalar ./(GA75*areascalar) * perareascalar;
GNMR = GNM.*co2scalar ./(GNA*areascalar) * perareascalar;
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,1) = cellstr(opporder(oo));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,2) = num2cell(round(GTA(:,2)));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,3) = num2cell(round(GTC(:,2)));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,4) = cellstr(strcat(num2str(round(GTM(:,2)),'%i')," [",...
    num2str(round(GTM(:,3)),'%i')," - ",num2str(round(GTM(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,5) = cellstr(strcat(num2str(round(GTMR(:,2)),'%i')," [",...
    num2str(round(GTMR(:,3)),'%i')," - ",num2str(round(GTMR(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,6) = cellstr(strcat(num2str(round(GTMO(:,2)),'%i')," [",...
    num2str(round(GTMO(:,3)),'%i')," - ",num2str(round(GTMO(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,7) = cellstr(strcat(num2str(round(GRPA(:,2)),'%i')," [",...
    num2str(round(GRPA(:,3)),'%i')," - ",num2str(round(GRPA(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,8) = cellstr(strcat(num2str(round(GPM(:,2)),'%i')," [",...
    num2str(round(GPM(:,3)),'%i')," - ",num2str(round(GPM(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,9) = cellstr(strcat(num2str(round(GPMR(:,2)),'%i')," [",...
    num2str(round(GPMR(:,3)),'%i')," - ",num2str(round(GPMR(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,10) = cellstr(strcat(num2str(round(GRA75(:,2)),'%i')," [",...
    num2str(round(GRA75(:,3)),'%i')," - ",num2str(round(GRA75(:,1)),'%i'),"]"));
T(nbiomes*(no+1)+3+1:nbiomes*(no+1)+3+no,11) = cellstr(strcat(num2str(round(GRNA(:,2)),'%i')," [",...
    num2str(round(GRNA(:,3)),'%i')," - ",num2str(round(GRNA(:,1)),'%i'),"]"));


% for bb = 1 : nbiomes
%     T((bb-1)*(no+1)+3,1) = cellstr(biomenames(bb));
%     for oo = 1 : no
%         T((bb-1)*(no+1)+3+oo,1) = cellstr(opporder(oo));
%         T((bb-1)*(no+1)+3+oo,2) = num2cell(round(TA(bb,oo,2),2));
%         T((bb-1)*(no+1)+3+oo,3) = num2cell(round(TC(bb,oo,2),2));
%         T((bb-1)*(no+1)+3+oo,4) = cellstr(strcat(num2str(round(TM(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(TM(bb,oo,3),2),'%.2f')," - ",num2str(round(TM(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,5) = cellstr(strcat(num2str(round(TMR(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(TMR(bb,oo,3),2),'%.2f')," - ",num2str(round(TMR(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,6) = cellstr(strcat(num2str(round(TMO(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(TMO(bb,oo,3),2),'%.2f')," - ",num2str(round(TMO(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,7) = cellstr(strcat(num2str(round(RPA(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(RPA(bb,oo,3),2),'%.2f')," - ",num2str(round(RPA(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,8) = cellstr(strcat(num2str(round(PM(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(PM(bb,oo,3),2),'%.2f')," - ",num2str(round(PM(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,9) = cellstr(strcat(num2str(round(PMR(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(PMR(bb,oo,3),2),'%.2f')," - ",num2str(round(PMR(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,10) = cellstr(strcat(num2str(round(RA75(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(RA75(bb,oo,3),2),'%.2f')," - ",num2str(round(RA75(bb,oo,1),2),'%.2f'),"]"));
%         T((bb-1)*(no+1)+3+oo,11) = cellstr(strcat(num2str(round(RNA(bb,oo,2),2),'%.2f')," [",...
%             num2str(round(RNA(bb,oo,3),2),'%.2f')," - ",num2str(round(RNA(bb,oo,1),2),'%.2f'),"]"));
%     end
% end

TT = cell2table(T);
writetable(TT,xcelfile,'Sheet','combined2')

% By Opportunity
% **************
oclm = ["Biomes","Area","No-albedo carbon",...
    "Total NCI (carbon+albedo)","Total NCI low","Total NCI (high)","uncertainty",...
    "NCI average rate","NCI rate low","NCI rate high",...
    "Mean Albedo Offset","Mean albedo low","Mean Albedo high",...
    "<50%AO Areas","<50%AO %Areas","% area low","% area high","uncertainty",...
    "<50%AO NCI","NCI low","NCI high","uncertainty","<50%AO NCI average rate",...
    ">75%AO Areas",">75%AO %Areas",">75%AO NCI",">75%AO NCI average rate",...
    ">100%AO Areas",">100%AO %Areas",">100%AO NCI",">100%AO NCI average rate"];

oclu = ["","[Mha]","[Pg CO2e]",...
    "[Pg CO2e]","[Pg CO2e]","[Pg CO2e]","[%]","[Mg/ha]","[Mg/ha]","[Mg/ha]",...
    "[%]","[%]","[%]","[Mha]","[%]","[%]","[%]","[%]",...
    "[Pg CO2e]","[Pg CO2e]","[Pg CO2e]","[%]","[Mg/ha]",...
    "[Mha]","[%]","[Pg CO2e]","[Mg/ha]","[Mha]","[%]","[Pg CO2e]","[Mg/ha]"];


for oo = 1 : no
    opp = opporder(oo);
    worldvalues = zeros(nbiomes,numel(oclm)-1);
    worldvalues(:,1) = TA(:,oo,2);
    worldvalues(:,2) = TC(:,oo,2);
    worldvalues(:,3) = TM(:,oo,2);
    worldvalues(:,4) = TM(:,oo,3);
    worldvalues(:,5) = TM(:,oo,1);
    worldvalues(:,6) = ((TM(:,oo,2)-TM(:,oo,3))./TM(:,oo,2)+(TM(:,oo,1)-TM(:,oo,2))./TM(:,oo,2))./2 .* 100;
    worldvalues(:,7) = TMR(:,oo,2);
    worldvalues(:,8) = TMR(:,oo,3);
    worldvalues(:,9) = TMR(:,oo,1);
    worldvalues(:,10) = TMO(:,oo,2);
    worldvalues(:,11) = TMO(:,oo,3);
    worldvalues(:,12) = TMO(:,oo,1);
    worldvalues(:,13) = PA(:,oo,2);
    worldvalues(:,14) = RPA(:,oo,2);
    worldvalues(:,15) = RPA(:,oo,3);
    worldvalues(:,16) = RPA(:,oo,1);
    worldvalues(:,17) = ((PA(:,oo,2)-PA(:,oo,3))./PA(:,oo,2)+(PA(:,oo,1)-PA(:,oo,2))./PA(:,oo,2))./2 .* 100;
    worldvalues(:,18) = PM(:,oo,2);
    worldvalues(:,19) = PM(:,oo,3);
    worldvalues(:,20) = PM(:,oo,1);
    worldvalues(:,21) = ((PM(:,oo,2)-PM(:,oo,3))./PM(:,oo,2)+(PM(:,oo,1)-PM(:,oo,2))./PM(:,oo,2))./2 .* 100;
    worldvalues(:,22) = PMR(:,oo,2);
    worldvalues(:,23) = A75(:,oo,2);
    worldvalues(:,24) = RA75(:,oo,2);
    worldvalues(:,25) = M75(:,oo,2);
    worldvalues(:,26) = M75R(:,oo,2);
    worldvalues(:,27) = NA(:,oo,2);
    worldvalues(:,28) = RNA(:,oo,2);
    worldvalues(:,29) = NM(:,oo,2);
    worldvalues(:,30) = NMR(:,oo,2);

    T = cell(nbiomes+2,numel(oclm));
    T(1,:) = cellstr(oclm);
    T(2,:) = cellstr(oclu);
    T(3:nbiomes+2,1) = cellstr(biomenames(semigeographicorder));
    T(3:nbiomes+2,2:numel(oclm)) = num2cell(worldvalues);
    
    TT = cell2table(T);
    writetable(TT,xcelfile,'Sheet',opp)    
end




