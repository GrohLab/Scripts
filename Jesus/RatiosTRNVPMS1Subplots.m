AllL200S1=ClustersL200OnsetS1(:);AllL200TRN=ClustersL200OnsetTRN(:); AllL200VPM=ClustersL200OnsetVPM(:);
AllControlOnS1=ClustersControlOnsetS1(:);AllControlOnTRN=ClustersControlOnsetTRN(:); AllControlOnVPM=ClustersControlOnsetVPM(:);

distL200OnvsControlOnS1=(hist(AllL200S1,1400)./300)./(hist(AllControlOnS1,1400)./150);
distL200OnvsControlOnS1=(distL200OnvsControlOnS1-1)*100;

distL200OnvsControlOnTRN=(hist(AllL200TRN,1400)./300)./(hist(AllControlOnTRN,1400)./150);
 distL200OnvsControlOnTRN=(distL200OnvsControlOnTRN-1)*100

distL200OnvsControlOnVPM=(hist(AllL200VPM,1400)./300)./(hist(AllControlOnVPM,1400)./150);(distL200OnvsControlOnVPM-1)*100
 distL200OnvsControlOnVPM=(distL200OnvsControlOnVPM-1)*100
figure;
subplot(4,1,1);
barPositions=1:size(distL200OnvsControlOnS1,2);
pFl=distL200OnvsControlOnS1>=0; NFl=distL200OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL200OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Potentiation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{1ms}'); ylabel('Ratio %'); 
title('Ratio AfterInduction/Control respone S1')
xticklabels({'-0.2', '0','0.2','0.4','0.6','0.8','1','1.2',}); hold off

subplot(4,1,2);
barPositions=1:size(distL200OnvsControlOnTRN,2);
pFl=distL200OnvsControlOnTRN>=0; NFl=distL200OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL200OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Potentiation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{1ms}'); ylabel('Ratio %'); 
title('Ratio AferInduction/Control respone TRN')
xticklabels({'-0.2', '0','0.2','0.4','0.6','0.8','1','1.2',}); hold off

subplot(4,1,3);
barPositions=1:size(distL200OnvsControlOnVPM,2);
pFl=distL200OnvsControlOnVPM>=0; NFl=distL200OnvsControlOnVPM<0;
bar(barPositions(pFl),distL200OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Potentiation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{1ms}'); ylabel('Ratio %'); 
title('Ratio AfterInduction/Control respone VPM')
xticklabels({'-0.2', '0','0.2','0.4','0.6','0.8','1','1.2',}); hold off
saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/ratiosCircuitS1_TRN_VPM']);
subplot(4,1,4);
plot(TimePsth,histcounts(ClustersControlOnsetVPM,1400));xlim([-200 1200])

%% only on
distL200OnvsControlOnS1=hist(AllL200S1,800)-hist(AllControlOnS1,800);
distL200OnvsControlOnTRN=hist(AllL200TRN,800)-hist(AllControlOnTRN,800);
distL200OnvsControlOnVPM=hist(AllL200VPM,800)-hist(AllControlOnVPM,800);
figure;
subplot(3,1,1);
barPositions=1:size(distL200OnvsControlOnS1,2);
pFl=distL200OnvsControlOnS1>=0; NFl=distL200OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL200OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]); hold off

subplot(3,1,2);
barPositions=1:size(distL200OnvsControlOnTRN,2);
pFl=distL200OnvsControlOnTRN>=0; NFl=distL200OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL200OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]);hold off

subplot(3,1,3);
barPositions=1:size(distL200OnvsControlOnVPM,2);
pFl=distL200OnvsControlOnVPM>=0; NFl=distL200OnvsControlOnVPM<0;
bar(barPositions(pFl),distL200OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'});xlim([450 600]);
%% off
distL200OnvsControlOnS1=hist(AllL200S1,800)-hist(AllControlOnS1,800);
distL200OnvsControlOnTRN=hist(AllL200TRN,800)-hist(AllControlOnTRN,800);
distL200OnvsControlOnVPM=hist(AllL200VPM,800)-hist(AllControlOnVPM,800);
figure;
subplot(3,1,1);
barPositions=1:size(distL200OnvsControlOnS1,2);
pFl=distL200OnvsControlOnS1>=0; NFl=distL200OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL200OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 OFF respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([650 800]); hold off

subplot(3,1,2);
barPositions=1:size(distL200OnvsControlOnTRN,2);
pFl=distL200OnvsControlOnTRN>=0; NFl=distL200OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL200OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 OFF respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([650 800]);hold off

subplot(3,1,3);
barPositions=1:size(distL200OnvsControlOnVPM,2);
pFl=distL200OnvsControlOnVPM>=0; NFl=distL200OnvsControlOnVPM<0;
bar(barPositions(pFl),distL200OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 OFF respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'});xlim([650 800]);