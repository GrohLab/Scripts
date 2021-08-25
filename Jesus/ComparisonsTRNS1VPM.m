%% control
ClustersControlOnsetVPM1=ClustersControlOnsetVPM(ClustersControlOnsetVPM>0.002 & ClustersControlOnsetVPM>0.05);

AllControlOnTRN=(ClustersControlOnsetTRN(:)*1000)+3;
AllControlOnVPM=(ClustersControlOnsetVPM(:)*1000)+3;
AllControlOnS1=(ClustersControlOnsetS1(:)*1000)+3;
%on an off
figure;
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllControlOnTRN,1400)))
hold on
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllControlOnVPM,1400)))
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllControlOnS1,1400)))
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Z-Score");
title("Circuit whisker Response Control"); 
%%onset
figure;
plot([-30:0.5:39.5],zscore(hist(AllControlOnTRN(AllControlOnTRN<40 & AllControlOnTRN>-30),140)));
hold on;
plot([-30:0.5:39.5],zscore(hist(AllControlOnVPM(AllControlOnVPM<40 & AllControlOnVPM>-30),140)));
plot([-30:0.5:39.5],zscore(hist(AllControlOnS1(AllControlOnS1<40 & AllControlOnS1>-30),140)));
title("Circuit Whisker Onset Response Control"); 
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Z-Score");

%figure all
figure;
histogram(AllControlOnTRN(AllControlOnTRN<80 & AllControlOnTRN>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);hold on
histogram(AllControlOnVPM(AllControlOnVPM<80 & AllControlOnVPM>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
histogram(AllControlOnS1(AllControlOnS1<80 & AllControlOnS1>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
title("Control");xlabel("time in ms");ylabel("probability");legend("TRN","VPM","S1")

figure; %on
histogram(AllControlOnTRN(AllControlOnTRN<40 & AllControlOnTRN>-5),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllControlOnVPM(AllControlOnVPM<40 & AllControlOnVPM>-5),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllControlOnS1(AllControlOnS1<40 & AllControlOnS1>-10),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Control On response");xlim([-5 40]);

figure; %Off
histogram(AllControlOnTRN(AllControlOnTRN<90 & AllControlOnTRN>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllControlOnVPM(AllControlOnVPM<90 & AllControlOnVPM>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllControlOnS1(AllControlOnS1<90 & AllControlOnS1>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Control Off response");

%%histogram(AllControlOnS1,400)))

%% L200

AllL200TRN=(ClustersL200OnsetTRN(:)*1000)+3;
AllL200VPM=(ClustersL200OnsetVPM(:)*1000)+3;
AllL200S1=(ClustersL200OnsetS1(:)*1000)+3;
%%on an off
figure;
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllL200TRN,1400)))
hold on
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllL200VPM,1400)))
plot([-0.2:0.001:1.2-0.001],zscore(hist(AllL200S1,1400)))
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Z-Score");
title("Circuit Response AfterInduction"); 
%%onset L200
figure;
plot([-100:1:49],zscore(hist(AllL200TRN(AllL200TRN<50 & AllL200TRN>-100),150)));
hold on;
plot([-100:1:49],zscore(hist(AllL200VPM(AllL200VPM<50 & AllL200VPM>-100),150)));
plot([-100:1:49],zscore(hist(AllL200S1(AllL200S1<50 & AllL200S1>-100),150)));
title("Circuit Whisker Onset Response + Laser Delay 200ms"); 
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Z-Score");
%figure all
figure;
histogram(AllL200TRN(AllL200TRN<80 & AllL200TRN>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);hold on
histogram(AllL200VPM(AllL200VPM<80 & AllL200VPM>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
histogram(AllL200S1(AllL200S1<80 & AllL200S1>-5),85,'Normalization','probability','DisplayStyle','stairs','LineWidth',2)
title("After Induction");xlabel("time in ms");ylabel("probability");legend("TRN","VPM","S1")

figure; %On
histogram(AllL200TRN(AllL200TRN<40 & AllL200TRN>-5),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllL200VPM(AllL200VPM<40 & AllL200VPM>-5),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllL200S1(AllL200S1<40 & AllL200S1>-5),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Circuit Whisker Onset Response After-Induction");xlim([-5 40]);

figure; %Off
histogram(AllL200TRN(AllL200TRN<80 & AllL200TRN>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllL200VPM(AllL200VPM<80 & AllL200VPM>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllL200S1(AllL200S1<80 & AllL200S1>45),90,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Circuit Whisker Offset Response + Laser Delay 200ms")%;xlim([30 150]);
%% ONL1
%%on an off
AllL1TRN=(ClustersL1OnsetTRN(:)*1000)+3;
AllL1VPM=(ClustersL1OnsetVPM(:)*1000)+3;
AllL1S1=(ClustersL1OnsetS1(:)*1000)+3;
figure;
plot([-0.25:0.001:0.149],zscore(hist(AllL1TRN,400)))
hold on
plot([-0.25:0.001:0.149],zscore(hist(AllL1VPM,400)))
plot([-0.25:0.001:0.149],zscore(hist(AllL1S1,400)))
legend("TRN","VPM","S1")
title("Circuit whisker Response + Laser Delay 1ms"); 

%plot onset response
figure;
plot([-100:1:49],zscore(hist(AllL1TRN(AllL1TRN<50 & AllL1TRN>-100),150)));
hold on;
plot([-100:1:49],zscore(hist(AllL1VPM(AllL1VPM<50 & AllL1VPM>-100),150)));
plot([-100:1:49],zscore(hist(AllL1S1(AllL1S1<50 & AllL1S1>-100),150)));
title("Circuit Onset Response"); 
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Z-Score");
title("Circuit Whisker Onset Response + Laser Delay 1ms"); 

figure; %On
histogram(AllL1TRN(AllL1TRN<50 & AllL1TRN>-80),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllL1VPM(AllL1VPM<50 & AllL1VPM>-80),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllL1S1(AllL1S1<50 & AllL1S1>-80),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Circuit Whisker Onset Response + Laser Delay -1ms");xlim([-70 50]);

figure; %Off
histogram(AllL1TRN(AllL1TRN<150 & AllL1TRN>30),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
hold on;
histogram(AllL1VPM(AllL1VPM<150 & AllL1VPM>30),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
histogram(AllL1S1(AllL1S1<150 & AllL1S1>30),150,'Normalization','probability','DisplayStyle','stairs','LineWidth',2);
legend("TRN","VPM","S1"); xlabel ("Time in ms"); ylabel("Probability")
title("Circuit Whisker Offset Response + Laser Delay -1ms");xlim([30 150]);
%% boxplot 

circuitbox=NaN(200000,3);
circuitbox(1:size(AllControlOnTRN(AllControlOnTRN<20 & AllControlOnTRN>3),1),1)= AllControlOnTRN(AllControlOnTRN<20 & AllControlOnTRN>3);
circuitbox(1:size(AllControlOnVPM(AllControlOnVPM<20 & AllControlOnVPM>3),1),2)=AllControlOnVPM(AllControlOnVPM<20 & AllControlOnVPM>3);
circuitbox(1:size(AllControlOnS1(AllControlOnS1<20 & AllControlOnS1>3),1),3)=AllControlOnS1(AllControlOnS1<20 & AllControlOnS1>3);

a= AllControlOnTRN(AllControlOnTRN<20 & AllControlOnTRN>3);%TRN
b= AllControlOnVPM(AllControlOnVPM<20 & AllControlOnVPM>3);%VPM
c= AllControlOnS1(AllControlOnS1<20 & AllControlOnS1>3);%S1
Onmarks=["TRN","VPM","S1"];
figure; subplot(1,2,1);
cdfplot(a); hold on; cdfplot(b);cdfplot(c); ylabel("probability");
legend("TRN","VPM","S1");subplot(1,2,2);
bar(mode(c))
hold on
bar(mode(a))
bar(mode(b))
legend("S1","TRN","VPM"); title("most frequent value (Mode)");ylabel("time(ms)");
%% TRN VPM and S1 IMAGESC diff and Zsscore
%diff imagesc S1
diffS1=hist(ClustersL200OnsetS1,400)-hist(ClustersControlOnsetS1,400);
diffL1S1=hist(ClustersL1OnsetS1,400)-hist(ClustersControlOnsetS1,400);
figure(10);
Ontimes=[-250:1:149];
yS1=1:size(diffS1,2);
imagesc(yS1,Ontimes,diffS1);colormap(jet); %L200-control
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-30 30]);
title("S1 diff Whiskers - whiskers+Laser delay -200ms");xlabel("Clusters");ylabel("time in ms")
figure(20); imagesc(yS1,Ontimes,diffL1S1);colormap(jet);%L1-control
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-10 10]);
title("S1 diff Whiskers - whiskers+Laser delay -1ms");xlabel("Clusters");ylabel("time in ms")
%TRN diff imagesc
diffTRN=hist(ClustersL200OnsetTRN,400)-hist(ClustersControlOnsetTRN,400);
diffL1TRN=hist(ClustersL1OnsetTRN,400)-hist(ClustersControlOnsetTRN,400);
figure(11);
Ontimes=[-250:1:149];
yTRN=1:size(diffTRN,2);
imagesc(yTRN,Ontimes,diffTRN);colormap(jet); %L200-control
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-10 10]);
title("TRN diff Whiskers - whiskers+Laser delay -200ms");xlabel("Clusters");ylabel("time in ms")
figure(21);imagesc(yTRN,Ontimes,diffL1TRN);colormap(jet);
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-35 35]);
title("TRN diff Whiskers - whiskers+Laser delay -1ms");xlabel("Clusters");ylabel("time in ms")
%VPM diff imagesc
diffVPM=hist(ClustersL200OnsetVPM,400)-hist(ClustersControlOnsetVPM,400);
diffL1VPM=hist(ClustersL1OnsetVPM,400)-hist(ClustersControlOnsetVPM,400);
figure(12); %L200-control
Ontimes=[-250:1:149];
yVPM=1:size(diffVPM,2);
imagesc(yVPM,Ontimes,diffVPM);colormap(jet);
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-10 10]);
title("VPM diff Whiskers - whiskers+Laser delay -200ms");xlabel("Clusters");ylabel("time in ms")
figure(22);%L1-Control
imagesc(yVPM,Ontimes,diffL1VPM);colormap(jet);
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-10 10]);
title("VPM diff Whiskers - whiskers+Laser delay -1ms");xlabel("Clusters");ylabel("time in ms")

%imagesc psth per clusters on response VPM
CircuitOnVPM = NaN(2000, size(ClustersControlOnsetVPM,2));
for i=1:size(ClustersControlOnsetVPM,2)
   value = ClustersControlOnsetVPM(:,i)<0.050 & ClustersControlOnsetVPM(:,i)>-0.070;
   spike=(ClustersControlOnsetVPM(value,i).*1000)+3;
    CircuitOnVPM(1:size(spike,1),i)=spike;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnVPM=hist(CircuitOnVPM,130);
yVPM=1:size(CircuitPsthOnVPM,2);
imagesc(yVPM,Ontimes,zscore(CircuitPsthOnVPM));colormap(jet);
colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-5 5]);
title("Control VPM");xlabel("Clusters");ylabel("time in ms")

%imagesc psth per clusters on response TRN
CircuitOnTRN = NaN(2000, size(ClustersControlOnsetTRN,2));
for i=1:size(ClustersControlOnsetTRN,2)
   value = ClustersControlOnsetTRN(:,i)<0.050 & ClustersControlOnsetTRN(:,i)>-0.070;
   spike=(ClustersControlOnsetTRN(value,i).*1000)+3;
    CircuitOnTRN(1:size(spike,1),i)=spike;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnTRN=hist(CircuitOnTRN,130);
yTRN=1:size(CircuitPsthOnTRN,2);
imagesc(yTRN,Ontimes,zscore(CircuitPsthOnTRN));colormap(jet);colorbar; ax=gca; ax.XAxisLocation='top';
clim=caxis; caxis([-4 4]); title("Control TRN");xlabel("Clusters");ylabel("time in ms")
%S1
CircuitOnS1 = NaN(2000, size(ClustersControlOnsetS1,2));
for i=1:size(ClustersControlOnsetS1,2)
   value = ClustersControlOnsetS1(:,i)<0.050 & ClustersControlOnsetS1(:,i)>-0.070;
   spikeS1=(ClustersControlOnsetS1(value,i).*1000)+3;
    CircuitOnS1(1:size(spikeS1,1),i)=spikeS1;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnS1=hist(CircuitOnS1,130);
yS1=1:size(CircuitPsthOnS1,2);
imagesc(yS1,Ontimes,zscore(CircuitPsthOnS1));colormap(jet);colorbar; ax=gca; ax.XAxisLocation='top';
clim=caxis; caxis([-5 5]); title("Control S1");xlabel("Clusters");ylabel("time in ms")

%% TRN VPM and S1 IMAGESC L200 Zsscore
%imagesc psth per clusters on response VPM
CircuitOnL200VPM = NaN(20000, size(ClustersL200OnsetVPM,2));
for i=1:size(ClustersL200OnsetVPM,2)
   value = ClustersL200OnsetVPM(:,i)<0.050 & ClustersL200OnsetVPM(:,i)>-0.070;
   spike=(ClustersL200OnsetVPM(value,i).*1000)+3;
    CircuitOnL200VPM(1:size(spike,1),i)=spike;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnL200VPM=hist(CircuitOnL200VPM,130);
yVPM=1:size(CircuitPsthOnL200VPM,2);
imagesc(yVPM,Ontimes,zscore(CircuitPsthOnL200VPM));
colormap(jet);colorbar; ax=gca; ax.XAxisLocation='top'; clim=caxis; caxis([-5 5]);
title("Whiskers + Laser delay 200ms VPM");xlabel("Clusters");ylabel("time in ms")

%imagesc psth per clusters on response TRN
CircuitOnL200TRN = NaN(20000, size(ClustersL200OnsetTRN,2));
for i=1:size(ClustersL200OnsetTRN,2)
   value = ClustersL200OnsetTRN(:,i)<0.050 & ClustersL200OnsetTRN(:,i)>-0.070;
   spikeL200TRN=(ClustersL200OnsetTRN(value,i).*1000)+3;
    CircuitOnL200TRN(1:size(spikeL200TRN,1),i)=spikeL200TRN;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnL200TRN=hist(CircuitOnL200TRN,130);
yTRN=1:size(CircuitPsthOnL200TRN,2);
imagesc(yTRN,Ontimes,zscore(CircuitPsthOnL200TRN));colormap(jet);colorbar; ax=gca; ax.XAxisLocation='top';
clim=caxis; caxis([-5 5]);
%S1
CircuitOnL200S1 = NaN(2000, size(ClustersL200OnsetS1,2));
for i=1:size(ClustersL200OnsetS1,2)
   value = ClustersL200OnsetS1(:,i)<0.050 & ClustersL200OnsetS1(:,i)>-0.070;
   spikeL200S1=(ClustersL200OnsetS1(value,i).*1000)+3;
    CircuitOnL200S1(1:size(spikeL200S1,1),i)=spikeL200S1;
end
figure;
Ontimes=[-70:1:49];
CircuitPsthOnL200S1=hist(CircuitOnL200S1,130);
yS1=1:size(CircuitPsthOnL200S1,2);
imagesc(yS1,Ontimes,zscore(CircuitPsthOnL200S1));colormap(jet);colorbar; ax=gca; ax.XAxisLocation='top';
clim=caxis; caxis([-5 5]);