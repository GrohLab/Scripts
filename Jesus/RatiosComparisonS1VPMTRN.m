%S1 L200
distL200OnvsControlOnS1=hist(AllL200S1,800)-hist(AllControlOnS1,800);
distL1OnvsControlOnS1=hist(AllL1S1,800)-hist(AllControlOnS1,800);
barPositions=1:size(distL200OnvsControlOnS1,2);
pFl=distL200OnvsControlOnS1>=0; NFl=distL200OnvsControlOnS1<0;
%L200-control s1
figure(1); bar(barPositions(pFl),distL200OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.emf']);
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.pdf']);
%L1
barPositions=1:size(distL1OnvsControlOnS1,2);
pFl=distL1OnvsControlOnS1>=0; NFl=distL1OnvsControlOnS1<0;
figure(2);
bar(barPositions(pFl),distL1OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL1OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser1 On response S1')
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.emf']);
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.pdf']);
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%% TRN
distL200OnvsControlOnTRN=hist(AllL200TRN,800)-hist(AllControlOnTRN,800);
distL1OnvsControlOnTRN=hist(AllL1TRN,800)-hist(AllControlOnTRN,800);
barPositions=1:size(distL200OnvsControlOnTRN,2);
pFl=distL200OnvsControlOnTRN>=0; NFl=distL200OnvsControlOnTRN<0; %%L200
figure(3);  bar(barPositions(pFl),distL200OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.emf']);
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.pdf']);
barPositions=1:size(distL1OnvsControlOnTRN,2);
pFl=distL1OnvsControlOnTRN>=0; NFl=distL1OnvsControlOnTRN<0;
figure(4);
bar(barPositions(pFl),distL1OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL1OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control-Laser1 On response TRN')
%saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.emf']);
%saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.pdf']);
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%% VPM
distL200OnvsControlOnVPM=hist(AllL200VPM,800)-hist(AllControlOnVPM,800);
distL1OnvsControlOnVPM=hist(AllL1VPM,800)-hist(AllControlOnVPM,800);
barPositions=1:size(distL200OnvsControlOnVPM,2);
pFl=distL200OnvsControlOnVPM>=0; NFl=distL200OnvsControlOnVPM<0;
figure(5);
bar(barPositions(pFl),distL200OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL200OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Ratio Control-Laser200 On respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.emf']);
%saveas(figure(1),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL200On.pdf']);
barPositions=1:size(distL1OnvsControlOnVPM,2);
pFl=distL1OnvsControlOnVPM>=0; NFl=distL1OnvsControlOnVPM<0;
figure(6);
bar(barPositions(pFl),distL1OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; bar(barPositions(NFl),distL1OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('different Number of Spikes'); 
title('Ratio Control-Laser1 On response VPM')
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.emf']);
saveas(figure(5),[pwd '/New figures Jesus/RatiosOn/ControlOnvsL1On.pdf']);
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
