% Ratio between Laser and Control
mkdir('New figures Jesus/RatiosOn')
mkdir('New figures Jesus/RatiosOff')
distL200OnvsControlOn=hist(allL200On,60)-hist(allControlOn,60);
distL200OffvsControlOff=hist(allL200Off,60)-hist(allControlOff,60);

distL100OnvsControlOn=hist(allL100On,60)-hist(allControlOn,60);
distL100OffvsControlOff=hist(allL100Off,60)-hist(allControlOff,60);

distL50OnvsControlOn=hist(allL50On,60)-hist(allControlOn,60);
distL50OffvsControlOff=hist(allL50Off,60)-hist(allControlOff,60);

distL10OnvsControlOn=hist(allL10On,60)-hist(allControlOn,60);
distL10OffvsControlOff=hist(allL10Off,60)-hist(allControlOff,60);

distL1OnvsControlOn=hist(allL1On,60)-hist(allControlOn,60);
distL1OffvsControlOff=hist(allL1Off,60)-hist(allControlOff,60);


barPositions=1:size(distL200OnvsControlOn,2);
pFl=distL200OnvsControlOn>=0; NFl=distL200OnvsControlOn<0;
figure(1); bar(barPositions(pFl),distL200OnvsControlOn(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOn(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser200 On respone')
saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.emf']);
saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.pdf']);

barPositions=1:size(distL100OnvsControlOn,2);
pFl=distL100OnvsControlOn>=0; NFl=distL100OnvsControlOn<0;
figure(2); bar(barPositions(pFl),distL100OnvsControlOn(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL100OnvsControlOn(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser100 On response')
saveas(figure(2),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL100On.emf']);
saveas(figure(2),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL100On.pdf']);

barPositions=1:size(distL50OnvsControlOn,2);
pFl=distL50OnvsControlOn>=0; NFl=distL50OnvsControlOn<0;
figure(3); 
bar(barPositions(pFl),distL50OnvsControlOn(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL50OnvsControlOn(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser50 On response')
saveas(figure(3),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL50On.emf']);
saveas(figure(3),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL50On.pdf']);

barPositions=1:size(distL10OnvsControlOn,2);
pFl=distL10OnvsControlOn>=0; NFl=distL10OnvsControlOn<0;
figure(4);
bar(barPositions(pFl),distL10OnvsControlOn(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL10OnvsControlOn(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser10 On response')
saveas(figure(4),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL10On.emf']);
saveas(figure(4),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL10On.pdf']);

barPositions=1:size(distL1OnvsControlOn,2);
pFl=distL1OnvsControlOn>=0; NFl=distL1OnvsControlOn<0;
figure(5);
bar(barPositions(pFl),distL1OnvsControlOn(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL1OnvsControlOn(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser1 On response')
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.emf']);
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.pdf']);
%off response ratio

barPositions=1:size(distL200OffvsControlOff,2);
pFl=distL200OffvsControlOff>=0; NFl=distL200OffvsControlOff<0;
figure(6); 
bar(barPositions(pFl),distL200OffvsControlOff(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL200OffvsControlOff(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser200 OFF response')
saveas(figure(6),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL200On.emf']);
saveas(figure(6),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL200On.pdf']);

barPositions=1:size(distL100OffvsControlOff,2);
pFl=distL100OffvsControlOff>=0; NFl=distL100OffvsControlOff<0;
figure(7); bar(barPositions(pFl),distL100OffvsControlOff(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL100OffvsControlOff(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser100 OFF response')
saveas(figure(7),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL100On.emf']);
saveas(figure(7),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL100On.pdf']);


barPositions=1:size(distL50OffvsControlOff,2);
pFl=distL50OffvsControlOff>=0; NFl=distL50OffvsControlOff<0;
figure(8); bar(barPositions(pFl),distL50OffvsControlOff(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL50OffvsControlOff(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser50 OFF response')
saveas(figure(8),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL50On.emf']);
saveas(figure(8),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL50On.pdf']);

barPositions=1:size(distL10OffvsControlOff,2);
pFl=distL10OffvsControlOff>=0; NFl=distL10OffvsControlOff<0;
figure(9); bar(barPositions(pFl),distL10OffvsControlOff(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL10OffvsControlOff(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser10 OFF response')

saveas(figure(9),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL10On.emf']);
saveas(figure(9),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL10On.pdf']);

barPositions=1:size(distL1OffvsControlOff,2);
pFl=distL1OffvsControlOff>=0; NFl=distL1OffvsControlOff<0;
figure(10); bar(barPositions(pFl),distL1OffvsControlOff(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL1OffvsControlOff(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control/Laser1 OFF response')
saveas(figure(10),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL1On.emf']);
saveas(figure(10),[pwd '/New figures Jesus/RatiosOff/ControlOnvsL1On.pdf']);

%% %Ratio On vs Off same condition





%% 
 
% Kurtosis On por clusters
figure(61);
kurtLControlOn =kurtosis(RspControlOn);
kurtL200On=kurtosis(RspL200On);
x=ones(size(kurtLControlOn)); x1=ones(size(kurtL200On))*2
plot([x;x1],[kurtLControlOn;kurtL200On],'-*');
xlim([0 3]);
title('Kurtosis ON response');
kurtallon=[kurtLControlOn' kurtL200On'];
boxplot(kurtallon,onLabels)
xlim([0 3])
%kurtosis off por clusters
figure(62);
kurtLControlOff=kurtosis(RspControlOff);
kurtL200Off=kurtosis(RspL200Off);
plot([x;x1],[kurtLControlOff;kurtL200Off],'-*');
xlim([0 3]);
title('Kurtosis Off response');
%ratios ;
figure
ratioOnOffControl=hist(allControlOn,60)-hist(allControlOff,60);
ratioOnOffL200=hist(allL200On,60)-hist(allL200Off,60);
bar(ratioOnOffControl)
bar(ratioONOffL200)
%%
figure;
skewness(allL1On);
figure;
skewness(allL1Off);
figure;
kurtosis(allL200On);
figure;
bar(kurtosis(allL200Off));

%%
figure;
subplot(1,2,1)
namesConditions={'Control' 'L200' '100' 'L50' 'L10' 'L1'};
DistMMControlOn=nanmean(allControlOn)-nanmedian(allControlOn);
DistMML200On=nanmean(allL200On)-nanmedian(allL200On);
DistMML100On=nanmean(allL100On)-nanmedian(allL100On);
DistMML50On=nanmean(allL50On)-nanmedian(allL50On);
DistMML10On=nanmean(allL10On)-nanmedian(allL10On);
DistMML1On=nanmean(allL1On)-nanmedian(allL1On);
DistMMOn=[DistMMControlOn DistMML200On DistMML100On DistMML50On DistMML10On DistMML1On];
bar(DistMMOn);set(gca,'XTickLabel',namesConditions)
ylabel('time _{(ms)}')
title ('Mean-Median Distance On')
subplot(1,2,2)
DistMMControlOff=nanmean(allControlOff)-nanmedian(allControlOff);
DistMML200Off=nanmean(allL200Off)-nanmedian(allL200Off);
DistMML100Off=nanmean(allL100Off)-nanmedian(allL100Off);
DistMML50Off=nanmean(allL50Off)-nanmedian(allL50Off);
DistMML10Off=nanmean(allL10Off)-nanmedian(allL10Off);
DistMML1Off=nanmean(allL1Off)-nanmedian(allL1Off);
DistMMOff=[DistMMControlOff DistMML200Off DistMML100Off DistMML50Off DistMML10Off DistMML1Off];
namesConditions={'Control' 'L200' '100' 'L50' 'L10' 'L1'};
bar(DistMMOff);set(gca,'XTickLabel',namesConditions);
ylabel('time _{(ms)}')
title ('Mean-Median Distance Off')

%% %Entropy
figure(100)
namesConditions={'Control','L200','L100','L50','L10','L1'}
EnControlOn=getEntropyFromPDF(hist(allControlOn,60));
EnL200On=getEntropyFromPDF(hist(allL200On,60));
EnL100On=getEntropyFromPDF(hist(allL100On,60));
EnL50On=getEntropyFromPDF(hist(allL50On,60));
EnL10On=getEntropyFromPDF(hist(allL10On,60));
EnL1On=getEntropyFromPDF(hist(allL1On,60));
EnOn=[EnControlOn EnL200On EnL100On EnL50On EnL10On EnL1On];
bar(EnOn); set(gca,'XtickLabel',namesConditions);
ylabel('Entropy')
ylim([0.7 1])
title('Entropy On response')
% OFF entropy
figure(101)
EnControlOff=getEntropyFromPDF(hist(allControlOff,60));
EnL200Off=getEntropyFromPDF(hist(allL200Off,60));
EnL100Off=getEntropyFromPDF(hist(allL100Off,60));
EnL50Off=getEntropyFromPDF(hist(allL50Off,60));
EnL10Off=getEntropyFromPDF(hist(allL10Off,60));
EnL1Off=getEntropyFromPDF(hist(allL1Off,60));
EnOff=[EnControlOff EnL200Off EnL100Off EnL50Off EnL10Off EnL1Off];
bar(EnOff); set(gca,'XtickLabel',namesConditions);
ylabel('Entropy')
ylim([0.7 1])
title ('Entropy off response')
saveas(figure(100),[pwd '/New figures Jesus/all/EntropyOn.emf']);
saveas(figure(101),[pwd '/New figures Jesus/all/EntropyOn.emf']);
saveas(figure(100),[pwd '/New figures Jesus/all/EntropyOff.pdf'])
saveas(figure(101),[pwd '/New figures Jesus/all/EntropyOff.pdf'])


