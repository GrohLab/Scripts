distL1OnvsControlOnS1=hist(AllL1S1,800)-hist(AllControlOnS1,800);
distL1OnvsControlOnTRN=hist(AllL1TRN,800)-hist(AllControlOnTRN,800);
distL1OnvsControlOnVPM=hist(AllL1VPM,800)-hist(AllControlOnVPM,800);

figure;
subplot(3,1,1);
barPositions=1:size(distL1OnvsControlOnS1,2);
pFl=distL1OnvsControlOnS1>=0; NFl=distL1OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL1OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('difference Whissker Control-Whisker&Laser delay -1ms On respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); hold off

subplot(3,1,2);
barPositions=1:size(distL1OnvsControlOnTRN,2);
pFl=distL1OnvsControlOnTRN>=0; NFl=distL1OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL1OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control-Whisker&Laser delay -1ms On respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); hold off

subplot(3,1,3);
barPositions=1:size(distL1OnvsControlOnVPM,2);
pFl=distL1OnvsControlOnVPM>=0; NFl=distL1OnvsControlOnVPM<0;
bar(barPositions(pFl),distL1OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('difference Whisker Control- Whisker&Laser delay -1ms On respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'})
%% only on
distL1OnvsControlOnS1=hist(AllL1S1,800)-hist(AllControlOnS1,800);
distL1OnvsControlOnTRN=hist(AllL1TRN,800)-hist(AllControlOnTRN,800);
distL1OnvsControlOnVPM=hist(AllL1VPM,800)-hist(AllControlOnVPM,800);

figure;
subplot(3,1,1);
barPositions=1:size(distL1OnvsControlOnS1,2);
pFl=distL1OnvsControlOnS1>=0; NFl=distL1OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL1OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('difference Whissker Control-Whisker&Laser delay -1ms On respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]); hold off

subplot(3,1,2);
barPositions=1:size(distL1OnvsControlOnTRN,2);
pFl=distL1OnvsControlOnTRN>=0; NFl=distL1OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL1OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control-Whisker&Laser delay -1ms On respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]); hold off

subplot(3,1,3);
barPositions=1:size(distL1OnvsControlOnVPM,2);
pFl=distL1OnvsControlOnVPM>=0; NFl=distL1OnvsControlOnVPM<0;
bar(barPositions(pFl),distL1OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('difference Whisker Control- Whisker&Laser delay -1ms On respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'});xlim([450 600]);
%% only on

figure;
subplot(3,1,1);
barPositions=1:size(distL1OnvsControlOnS1,2);
pFl=distL1OnvsControlOnS1>=0; NFl=distL1OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL1OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms ON respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]); hold off

subplot(3,1,2);
barPositions=1:size(distL1OnvsControlOnTRN,2);
pFl=distL1OnvsControlOnTRN>=0; NFl=distL1OnvsControlOnTRN<0; %%L200
%figure(3);
bar(barPositions(pFl),distL1OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms ON respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([450 600]);hold off

subplot(3,1,3);
barPositions=1:size(distL1OnvsControlOnVPM,2);
pFl=distL1OnvsControlOnVPM>=0; NFl=distL1OnvsControlOnVPM<0;
bar(barPositions(pFl),distL1OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('difference Number of Spikes'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms ON respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'});xlim([450 600]);
%% off

figure;
subplot(3,1,1);
barPositions=1:size(distL1OnvsControlOnS1,2);
pFl=distL1OnvsControlOnS1>=0; NFl=distL1OnvsControlOnS1<0;
%L200-control s1
bar(barPositions(pFl),distL1OnvsControlOnS1(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnS1(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms OFF respone S1')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([650 800]); hold off

subplot(3,1,2);
barPositions=1:size(distL1OnvsControlOnTRN,2);
pFl=distL1OnvsControlOnTRN>=0; NFl=distL1OnvsControlOnTRN<0; %%L1
%figure(3);
bar(barPositions(pFl),distL1OnvsControlOnTRN(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnTRN(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms OFF respone TRN')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'}); xlim([650 800]);hold off

subplot(3,1,3);
barPositions=1:size(distL1OnvsControlOnVPM,2);
pFl=distL1OnvsControlOnVPM>=0; NFl=distL1OnvsControlOnVPM<0;
bar(barPositions(pFl),distL1OnvsControlOnVPM(pFl),'FaceColor','g','DisplayName','Facilitation');
hold on; 
bar(barPositions(NFl),distL1OnvsControlOnVPM(NFl),'FaceColor','r','DisplayName','Depression');
legend
xticks('auto')
xlabel('bin_{0.5ms}'); ylabel('Counts'); 
title('Difference Whisker Control- Whisker&Laser delay -1ms OFF respone VPM')
xticklabels({'-250', '-200','150','-100','-50','0','50','100','150'});xlim([650 800]);