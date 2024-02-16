clear all

% ensemble from Joos et al. with 17th being ensem mean
% data extracted only to 200 years of the IRF
load ocean_uptake.mat
fremain = [1;irf(:,17);ones(199,1).*irf(200,17)];
yremain = 1:length(fremain);
GWP100 = sum(fremain(1:100))./100; % average CO2 concentration remaining

yremove = 0.01:1:200.01;
nremoveyears = 200;
nremainyears = 400;

% Carbon Accumulation Curve 1 - Chapman Richards
% example curve: https://www.sciencedirect.com/science/article/pii/S1674775518305109?via%3Dihub
Bmax=20; % kg C m-2 arbitrary magnitude chosen as representative value
kappa=.8; % keep between 0.5 and 1; 
mu = .05; % lower value delays rise, very strong sensitivity over range 0.05 to 0.1
lambda = 5; % higher value delays initial rise
text1 = {'Curve 1  k=.8, mu=.05, lambda=5'};
cumCremove1 = Bmax.*(1-kappa.*exp(-mu.*yremove)).^lambda;
Cremove1 = diff(cumCremove1);
RemoveRemain1=zeros(nremoveyears,nremainyears+nremoveyears);
for thisyremove = 1:nremoveyears
    for tremain=1:nremainyears
        temp = Cremove1(thisyremove).*fremain(tremain);
        RemoveRemain1(thisyremove,tremain+thisyremove) = temp;
    end
end
cumCremain1 = sum(RemoveRemain1,1);
% NOTE: cumCremain inaccurate after 200 years when our IRF curve stops
%       such that fremain is not decreasing after 200 years        


% Carbon Accumulation Curve 2 - Chapman Richards
Bmax=20;
kappa=.8; % keep between 0.5 and 1; 
mu = .03; % lower value delays rise, very strong sensitivity over range 0.05 to 0.1
lambda = 5; % higher value delays initial rise
text2 = {'Curve 2  k=.8, mu=.03, lambda=5'};
cumCremove2 = Bmax.*(1-kappa.*exp(-mu.*yremove)).^lambda;
Cremove2 = diff(cumCremove2);
RemoveRemain2=zeros(nremoveyears,nremainyears+nremoveyears);
for thisyremove = 1:nremoveyears
    for tremain=1:nremainyears
        temp = Cremove2(thisyremove).*fremain(tremain);
        RemoveRemain2(thisyremove,tremain+thisyremove) = temp;
    end
end
cumCremain2 = sum(RemoveRemain2,1);


% Carbon Accumulation Curve 3 - Chapman Richards
Bmax=20;
kappa=.8; % keep between 0.5 and 1; 
mu = .04; % lower value delays rise, very strong sensitivity over range 0.05 to 0.1
lambda = 6; % higher value delays initial rise
text3 = {'Curve 3  k=.8, mu=.04, lambda=6'};
cumCremove3 = Bmax.*(1-kappa.*exp(-mu.*yremove)).^lambda;
Cremove3 = diff(cumCremove3);
RemoveRemain3=zeros(nremoveyears,nremainyears+nremoveyears);
for thisyremove = 1:nremoveyears
    for tremain=1:nremainyears
        temp = Cremove3(thisyremove).*fremain(tremain);
        RemoveRemain3(thisyremove,tremain+thisyremove) = temp;
    end
end
cumCremain3 = sum(RemoveRemain3,1);


% average removal remaining from 25 to 100 years
avgwindow = 70:100;
fCremain1_avgwindow = mean(cumCremain1(avgwindow)./max(cumCremove1));
fCremain2_avgwindow = mean(cumCremain2(avgwindow)./max(cumCremove2));
fCremain3_avgwindow = mean(cumCremain3(avgwindow)./max(cumCremove3));


figure(1);clf
plot(yremain,fremain,'LineWidth',3)
ylabel('Fraction Remaining')
xlabel('Years after Pulse Removal')
fname=('FractionRemovalsRemaining.jpg'); print(fname,'-djpeg','-r300')

figure(2);clf
subplot(2,1,1)
plot(yremove(1:200),Cremove1,'b','LineWidth',3);hold on
plot(yremove(1:200),Cremove2,'c','LineWidth',3)
%plot(yremove(1:200),Cremove3,'k-','LineWidth',3)
ylabel({'Annual C Removal'; '[kg C m^-^2 y^-^1]'})
xlabel('year')
legend('Fast C Accumulation','Slow C Accumulation')
%legend([text1;text2;text3])
subplot(2,1,2)
plot(yremove,cumCremove1,'b-','LineWidth',3);hold on
plot(yremove,cumCremove2,'c-','LineWidth',3)
%plot(yremove,cumCremove3,'-k','LineWidth',3)
axis([0 200 0 20])
ylabel({'Cumulative Removal';'[kg C m^-^2]'})
xlabel('year')
fname=('CarbonRemovals.jpg'); print(fname,'-djpeg','-r300')

figure(3);clf
plot(RemoveRemain1(1,1:200),'m-','LineWidth',3);hold on
plot(RemoveRemain1(10,1:200),'y-','LineWidth',3)
plot(RemoveRemain1(30,1:200),'c-','LineWidth',3)
plot(RemoveRemain1(50,1:200),'g-','LineWidth',3)
plot(RemoveRemain1(100,1:200),'k-','LineWidth',3)
plot(yremove(1:200)+1,Cremove1,'b.')
legend('1 yr','10 yr','30 yr','50 yr','100 yr')
ylabel({'CO2e Mitigation from'; 'Year-Specific CO2 Removals Pulses';'[kg C m^-^2 y^-^1]'})
xlabel('year')
fname=('CarbonRemovalsPulses.jpg'); print(fname,'-djpeg','-r300')

figure(4);clf
plot(cumCremain1(1:nremainyears)./max(cumCremove1),'b-','LineWidth',3);hold on
plot(cumCremain2(1:nremainyears)./max(cumCremove2),'r-','LineWidth',3)
plot(cumCremain3(1:nremainyears)./max(cumCremove3),'k-','LineWidth',3)
plot([avgwindow(1) avgwindow(length(avgwindow))],[fCremain1_avgwindow fCremain1_avgwindow],'b:','LineWidth',3)
plot([avgwindow(1) avgwindow(length(avgwindow))],[fCremain2_avgwindow fCremain2_avgwindow],'r:','LineWidth',3)
plot([avgwindow(1) avgwindow(length(avgwindow))],[fCremain3_avgwindow fCremain3_avgwindow],'k:','LineWidth',3)
ylabel({'Cumulative Carbon-only Mitigation'; '[Fraction of Max Removal]'})
xlabel('year')
fname=('CMitigationRemaining.jpg'); print(fname,'-djpeg','-r300')


% ------------------------------------
%  CASE 1 Albedo = 100% of Carbon
% ------------------------------------

albedodeduction = 1; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove1;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

oceanlandrelease1 = cumCremove1(1:nremoveyears) - cumCremain1(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);


figure(5);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove1,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease1(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove1)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land'; '[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 1: Fast C Accumulation, Albedo = 100% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','southeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case1.jpg';
print(fname,'-djpeg','-r300')


% ------------------------------------
%  CASE 2 Albedo = 50% of Carbon
% ------------------------------------

albedodeduction = 0.5; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove1;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');


oceanlandrelease1 = cumCremove1(1:nremoveyears) - cumCremain1(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);

figure(6);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove1,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease1(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove1)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land';'[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 2: Fast C Accumulation, Albedo = 50% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','southeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case2.jpg';
print(fname,'-djpeg','-r300')


% ------------------------------------
%  CASE 3 Albedo = 10% of Carbon
% ------------------------------------

albedodeduction = 0.1; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove1;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove1.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');


oceanlandrelease1 = cumCremove1(1:nremoveyears) - cumCremain1(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove1(1:nremoveyears)+oceanlandrelease1(1:nremoveyears);

figure(7);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove1,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease1(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove1)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land';'[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 3: Fast C Accumulation, Albedo = 10% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','northeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case3.jpg';
print(fname,'-djpeg','-r300')




% ------------------------------------
%  CASE 4 Albedo = 100% of Carbon
% ------------------------------------

albedodeduction = 1; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove2;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

oceanlandrelease2 = cumCremove2(1:nremoveyears) - cumCremain2(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);


figure(15);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove2,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease2(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove2)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land';'[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 4: Slow C Accumulation, Albedo = 100% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','southeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case4.jpg';
print(fname,'-djpeg','-r300')



% ------------------------------------
%  CASE 5 Albedo = 50% of Carbon
% ------------------------------------

albedodeduction = 0.5; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove2;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');


oceanlandrelease2 = cumCremove2(1:nremoveyears) - cumCremain2(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);

figure(16);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove2,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease2(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove2)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land';'[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 5: Slow C Accumulation, Albedo = 50% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','southeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case5.jpg';
print(fname,'-djpeg','-r300')


% ------------------------------------
%  CASE 6 Albedo = 10% of Carbon
% ------------------------------------

albedodeduction = 0.1; % X percent deduction in climate benefit from carbon caused by albedo effect
albedolead = 0.5;  % fraction of carbon mid-point time by which albedo leads

curve = cumCremove2;
thesequantiles = [0:.01:1];
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo1 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');

albedolead = 0.8;  
qarray = quantile(curve,thesequantiles);
for i = 1:length(qarray)
    temp = abs(qarray(i)-curve);
    [scrap, yarray(i)] = min(temp);
end
ynewarray = yarray.*albedolead;
temp = cumCremove2.*albedodeduction;
cumAlbedo2 = interp1(ynewarray,temp(yarray),[1:nremoveyears],'linear','extrap');


oceanlandrelease2 = cumCremove2(1:nremoveyears) - cumCremain2(1:nremoveyears);
netclimate1 = cumAlbedo1-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);
netclimate2 = cumAlbedo2-cumCremove2(1:nremoveyears)+oceanlandrelease2(1:nremoveyears);

figure(17);clf
subplot(2,1,1)
plot(1:nremoveyears,cumAlbedo1,'g-','LineWidth',3);hold on
plot(1:nremoveyears,cumAlbedo2,'g:','LineWidth',3);hold on
plot(yremove,-cumCremove2,'b-','LineWidth',3)
plot(1:nremoveyears,oceanlandrelease2(1:nremoveyears),'r-','LineWidth',3)
plot([0 201], [0 0],'k')
maxminy = max([max(cumAlbedo1) max(cumAlbedo2) max(cumCremove2)]);
axis([0 201 -maxminy maxminy])
ylabel({'Cumulative CO2e Flux'; 'Atmosphere to Land';'[kg C m^-^2]'})
legend('Albedo 50% Earlier','Albedo 80% Earlier','Carbon Removal','Ocean + Land Release','Location','southeast')
title('Case 6: Slow C Accumulation, Albedo = 10% Carbon')
subplot(2,1,2)
plot(1:nremoveyears,netclimate1,'k-','LineWidth',3);hold on
plot(1:nremoveyears,netclimate2,'k:','LineWidth',3);hold on
plot([0 201], [0 0],'k')
maxminy = max([max(abs(netclimate1)) max(abs(netclimate2))]);
axis([0 201 -maxminy maxminy])
legend('Albedo 50% Earlier','Albedo 80% Earlier','Location','northeast')
ylabel({'Net Cumulative CO2e Flux'; '[kg C m^-^2]'})
fname = 'netclimate_case6.jpg';
print(fname,'-djpeg','-r300')



