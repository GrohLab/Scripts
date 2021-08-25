%cdfplot(control)
figure
subplot(2,5,[1 2])
VPMControlON1st=AllControlOnVPM(AllControlOnVPM<=0.020 & AllControlOnVPM>=0.002);
TRNControlON1st=AllControlOnTRN(AllControlOnTRN<=0.020 & AllControlOnTRN>=0.002);
S1ControlON1st=AllControlOnS1(AllControlOnS1<=0.020 & AllControlOnS1>=0.002);
cdfplot(VPMControlON1st); hold on
cdfplot(TRNControlON1st);
cdfplot(S1ControlON1st);
title("First stimulus Control");xlabel("time in ms");ylabel("Cumulative probability");
legend("VPM","TRN","S1");
saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/CDF_COntrol_1st_S1_TRN_VPM']);

%cdfplot (afterinducton

subplot(2,5,[6 7])
VPMAfterON1st=AllL200VPM(AllL200VPM<=0.020 & AllL200VPM>=0.002);
TRNafterON1st=AllL200TRN(AllL200TRN<=0.020 & AllL200TRN>=0.002);
S1L2001st=AllL200S1(AllL200S1<=0.020 & AllL200S1>=0.002);
cdfplot(VPMAfterON1st);hold on;
cdfplot(TRNafterON1st);
cdfplot(S1L2001st);
title("First stimulus After-Induction");xlabel("time in ms");ylabel("Cumulative probability");
legend("VPM","TRN","S1");
%%Boxplot control
S1stControl= ClustersControlOnsetS1>=0 & ClustersControlOnsetS1<=0.040;
VPM1stControl=ClustersControlOnsetVPM>=0 & ClustersControlOnsetVPM<=0.040;
TRN1stControl=ClustersControlOnsetTRN>=0 & ClustersControlOnsetTRN<=0.040;
subplot(2,5,[4 5])
boxplot(CircuitMode,'Labels',{'VPM','TRN','S1Bf'},'Orientation','horizontal');xlim([-0.010 0.043]); title(" Spike Time Mode Control");
xlabel("time in s")
subplot(2,5,[9 10])
boxplot(CircuitModeAfter,'Labels',{'VPM','TRN','S1Bf'},'Orientation','horizontal');xlim([-0.010 0.043]); title(" Spike Time Mode Control");
xlabel("Time in ms")
saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/CDF_andBoxplotMode_control_and_after_1st_S1_TRN_VPM']);

CircuitMode= NaN(300,3); CircuitMode(1:size(VPMMode1st),1)=VPMMode1st;
CircuitMode(1:size(TRNMode1st),2)=TRNMode1st; CircuitMode(1:size(s1Mode1st),3)=s1Mode1st;
CircuitModeAfter= NaN(300,3); CircuitModeAfter(1:size(VPMMode1stAfter),1)=VPMMode1stAfter;
CircuitModeAfter(1:size(TRNMode1stAfter),2)=TRNMode1stAfter; CircuitModeAfter(1:size(s1Mode1stAfter),3)=s1Mode1stAfter;
figure;
boxplot([CircuitMode(:,1) CircuitModeAfter(:,1)],'Orientation','horizontal');xlim([0 0.02]);
figure;
boxplot([CircuitMode(:,2) CircuitModeAfter(:,2)],'Orientation','horizontal');xlim([0 0.02]);
yticklabels({"TRN-Control","TRN-After-Induction"})
figure;
boxplot([CircuitMode(:,3) CircuitModeAfter(:,3)],'Orientation','horizontal');xlim([0 0.02]);
yticklabels({"S1-Control","S1-After-Induction"})


%% sequential mapa control
%MOCK
ClustersControlOnset1stMOCK=nan(size(ClustersControlOnsetMOCK,1),size(ClustersControlOnsetMOCK,2));
id=find(ClustersControlOnsetMOCK>=0.003 & ClustersControlOnsetMOCK<=0.04);
for i=1:size(id);
ClustersControlOnset1stMOCK(id(i))=ClustersControlOnsetMOCK(id(i));
end
ClustersControlOnset1stMOCK=(ClustersControlOnset1stMOCK)+0.0025;
boxplot(mode(ClustersControlOnset1stMOCK))
%vpm

ClustersControlOnset1stVPM=nan(size(ClustersControlOnsetVPM,1),size(ClustersControlOnsetVPM,2));
id=find(ClustersControlOnsetVPM>=0.003 & ClustersControlOnsetVPM<=0.04);
for i=1:size(id);
ClustersControlOnset1stVPM(id(i))=ClustersControlOnsetVPM(id(i));
end
ClustersControlOnset1stVPM=(ClustersControlOnset1stVPM)+0.0025;
boxplot(mode(ClustersControlOnset1stVPM))
%TRN
ClustersControlOnset1stTRN=nan(size(ClustersControlOnsetTRN,1),size(ClustersControlOnsetTRN,2));
id=find(ClustersControlOnsetTRN>=0.003 & ClustersControlOnsetTRN<=0.04);
for i=1:size(id);
ClustersControlOnset1stTRN(id(i))=ClustersControlOnsetTRN(id(i));
end
ClustersControlOnset1stTRN=(ClustersControlOnset1stTRN)+0.0025;
boxplot(mode(ClustersControlOnset1stTRN))

%S1
ClustersControlOnset1stS1=nan(size(ClustersControlOnsetS1,1),size(ClustersControlOnsetS1,2));
id=find(ClustersControlOnsetS1>=0.005 & ClustersControlOnsetS1<=0.04);
for i=1:size(id);
ClustersControlOnset1stS1(id(i))=ClustersControlOnsetS1(id(i));
end
ClustersControlOnset1stS1=(ClustersControlOnset1stS1)+0.0025;
boxplot(mode(ClustersControlOnset1stS1))

VPMMode1st=mode(ClustersControlOnset1stVPM)';
s1Mode1st=mode(ClustersControlOnset1stS1)';
TRNMode1st=mode(ClustersControlOnset1stTRN)';
figure;
CircuitMode= NaN(300,3); CircuitMode(1:size(VPMMode1st),1)=VPMMode1st;
CircuitMode(1:size(TRNMode1st),2)=TRNMode1st; CircuitMode(1:size(s1Mode1st),3)=s1Mode1st;
boxplot(CircuitMode)
secuantialactivationVPMTRNS1map=cat(2,ClustersControlOnset1stVPM(:,4:54),ClustersControlOnset1stTRN(:,1:50),ClustersControlOnset1stS1(:,100:150));
figure(1);
subplot(13,1,[1 2 3 4]); ax2=subplot(13,1,[1 2 3 4])
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersControlOnset1stVPM(:),100)),'LineWidth',2);hold on; plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersControlOnset1stTRN(:),100)),'LineWidth',2)
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersControlOnset1stS1(:),100)),'LineWidth',2);hold off; xlim([-0.005 0.04])
title("Sequential activation Control"), xlabel("time in s");ylabel("zscore")
figure(1)
subplot(13,1,[6 7 8 9 10 11 12 13]); ax3=subplot(13,1,[6 7 8 9 10 11 12 13])
secuantialactivationVPMTRNS1Controlmap=cat(2,ClustersControlOnset1stVPM(:,4:54),ClustersControlOnset1stTRN(:,1:50),ClustersControlOnset1stS1(:,100:150));
NtotalCircuitControl=1:size(secuantialactivationVPMTRNS1Controlmap,2);
ZscoresequentialPsthControl=zscore(hist(secuantialactivationVPMTRNS1Controlmap,50))
imagesc([-0.010:0.001:0.05-0.001],NtotalCircuitControl,ZscoresequentialPsthControl'); axis xy
colormap(jet);caxis([0 5]);hold on; plot([0,size(NtotalCircuitControl,2)],[0.0;0.0],'LineStyle','--','LineWidth',2,'Color','w')
hold off; xlim([-0.005 0.04]); xlabel("time in s");ylabel("Circuit clusters")
hold on
plot([-0.005 0.04],[42;42],'LineStyle','--','LineWidth',1,'Color','w')
plot([-0.005 0.04],[91;91],'LineStyle','--','LineWidth',1,'Color','w')
hold off

linkaxes([ax2,ax3],'x'); colorbar; titlecolorbar=colorbar; titlecolorbar.Label.String= "zscore"
saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/Sequential_Activation1st_S1_TRN_VPM_Control']);


%% mapa afterinduction sequential
ClustersL200Onset1stMOCK=nan(size(ClustersL200OnsetMOCK,1),size(ClustersL200OnsetMOCK,2));
id=find(ClustersL200OnsetMOCK>=0.003 & ClustersL200OnsetMOCK<=0.04);
for i=1:size(id);
    ClustersL200Onset1stMOCK(id(i))=ClustersL200OnsetMOCK(id(i));
end
ClustersL200Onset1stMOCK=(ClustersL200Onset1stMOCK)+0.0025;
boxplot(mode(ClustersL200Onset1stMOCK))


ClustersL200Onset1stVPM=nan(size(ClustersL200OnsetVPM,1),size(ClustersL200OnsetVPM,2));
id=find(ClustersL200OnsetVPM>=0.003 & ClustersL200OnsetVPM<=0.04);
for i=1:size(id);
    ClustersL200Onset1stVPM(id(i))=ClustersL200OnsetVPM(id(i));
end
ClustersL200Onset1stVPM=(ClustersL200Onset1stVPM)+0.0025;
boxplot(mode(ClustersL200Onset1stVPM))

ClustersL200Onset1stTRN=nan(size(ClustersL200OnsetTRN,1),size(ClustersL200OnsetTRN,2));
id=find(ClustersL200OnsetTRN>=0 & ClustersL200OnsetTRN<=0.04);
for i=1:size(id);
ClustersL200Onset1stTRN(id(i))=ClustersL200OnsetTRN(id(i));
end
ClustersL200Onset1stTRN=(ClustersL200Onset1stTRN)+0.0025;
boxplot(mode(ClustersL200Onset1stTRN))

ClustersL200Onset1stS1=nan(size(ClustersL200OnsetS1,1),size(ClustersL200OnsetS1,2));
id=find(ClustersL200OnsetS1>=0.003 & ClustersL200OnsetS1<=0.04);
for i=1:size(id);
ClustersL200Onset1stS1(id(i))=ClustersL200OnsetS1(id(i));
end
ClustersL200Onset1stS1=(ClustersL200Onset1stS1)+0.0025;
boxplot(mode(ClustersL200Onset1stS1))

VPMMode1stAfter=mode(ClustersL200Onset1stVPM)';
s1Mode1stAfter=mode(ClustersL200Onset1stS1)';
TRNMode1stAfter=mode(ClustersL200Onset1stTRN)';

subplot(13,1,[1 2 3 4]); ax2=subplot(13,1,[1 2 3 4])
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stVPM(:),100)),'LineWidth',2);hold on; plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stTRN(:),100)),'LineWidth',2)
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stS1(:),100)),'LineWidth',2);hold off; xlim([-0.005 0.04])
title("Sequential activation After-Induction"), xlabel("time in s");ylabel("zscore")

subplot(13,1,[6 7 8 9 10 11 12 13]); ax3=subplot(13,1,[6 7 8 9 10 11 12 13])
secuantialactivationVPMTRNS1AFTERmap=cat(2,ClustersL200Onset1stVPM(:,4:54),ClustersL200Onset1stTRN(:,1:50),ClustersL200Onset1stS1(:,100:150));
NtotalCircuit=1:size(secuantialactivationVPMTRNS1AFTERmap,2);
ZscoresequentialPsth=zscore(hist(secuantialactivationVPMTRNS1AFTERmap,50))
imagesc([-0.017:0.0001:0.053-0.001],NtotalCircuit,ZscoresequentialPsth'); axis xy
colormap(jet);caxis([-1 8]);hold on; plot([0;size(NtotalCircuit,2)],[0.0;0.0],'LineStyle','--','LineWidth',2,'Color','w')
hold off; xlim([-0.005 0.04]); xlabel("time in s");ylabel("Circuit clusters")

linkaxes([ax2,ax3],'x'); colorbar; titlecolorbar=colorbar; 
saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/Sequential_Activation1st_S1_TRN_VPM_After']);
%% probability diff 
%VPM
MatrixVPMProbAfter1st=histcounts(ClustersL200Onset1stVPM,50,'Normalization','probability')
MatrixVPMProbControl1st=histcounts(ClustersControlOnset1stVPM,50,'Normalization','probability');
diffProb1stVPM=MatrixVPMProbAfter1st-MatrixVPMProbControl1st;
plot(diffProb1stVPM);
%TRN
MatrixTRNProbAfter1st=histcounts(ClustersL200Onset1stTRN,50,'Normalization','probability')
MatrixTRNProbControl1st=histcounts(ClustersControlOnset1stTRN,50,'Normalization','probability');
diffProb1stTRN=MatrixTRNProbAfter1st-MatrixTRNProbControl1st;
plot(diffProb1stTRN);
%S1
MatrixS1ProbAfter1st=histcounts(ClustersL200Onset1stS1,50,'Normalization','probability')
MatrixS1ProbControl1st=histcounts(ClustersControlOnset1stS1,50,'Normalization','probability');
diffProb1stS1=MatrixS1ProbAfter1st-MatrixS1ProbControl1st;
plot(diffProb1stS1);
% plot diff in probability per structura VPM TRN S1
plot([-0.007:0.001:0.043-0.001],diffProb1stVPM,'LineWidth',2)
hold on
plot([-0.007:0.001:0.043-0.001],diffProb1stTRN,'LineWidth',2)
plot([-0.007:0.001:0.043-0.001],diffProb1stS1,'LineWidth',2)
title("Probability Change")
 legend("VPM","TRN","S1"); hold off;
 
%mock get 1st trial control
ClustersControlOnset1stMOCK=nan(size(ClustersControlOnsetMOCK,1),size(ClustersControlOnsetMOCK,2));
id=find(ClustersControlOnsetMOCK>=-0.010 & ClustersControlOnsetMOCK<=0.04);
for i=1:size(id);
ClustersControlOnset1stMOCK(id(i))=ClustersControlOnsetMOCK(id(i));
end
ClustersControlOnset1stMOCK=(ClustersControlOnset1stMOCK)+0.003;
boxplot(mode(ClustersControlOnset1stMOCK))
%mock get 1st trial After
ClustersL200Onset1stMOCK=nan(size(ClustersL200OnsetMOCK,1),size(ClustersL200OnsetMOCK,2));
id=find(ClustersL200OnsetMOCK>=-0.010 & ClustersL200OnsetMOCK<=0.04);
for i=1:size(id);
ClustersL200Onset1stMOCK(id(i))=ClustersL200OnsetMOCK(id(i));
end
ClustersL200Onset1stMOCK=(ClustersL200Onset1stMOCK)+0.003;
boxplot(mode(ClustersL200Onset1stMOCK))

%Ratio
MatrixMOCKProbAfter1st=histcounts(ClustersL200Onset1stMOCK,50,'Normalization','probability')
MatrixMOCKProbControl1st=histcounts(ClustersControlOnset1stMOCK,50,'Normalization','probability');
diffProb1stMOCK=MatrixMOCKProbAfter1st-MatrixMOCKProbControl1st;
plot(diffProb1stMOCK);

ratioProb1stVPM=((MatrixVPMProbAfter1st./MatrixVPMProbControl1st)-1)*100;
ratioProb1stTRN=((MatrixTRNProbAfter1st./MatrixTRNProbControl1st)-1)*100;
ratioProb1stS1=((MatrixS1ProbAfter1st./MatrixS1ProbControl1st)-1)*100;
ratioProb1stMOCK=((MatrixMOCKProbAfter1st./MatrixMOCKProbControl1st)-1)*100;
figure
plot([-0.007:0.001:0.043-0.001],ratioProb1stVPM,'LineWidth',2); hold on;
plot([-0.007:0.001:0.043-0.001],ratioProb1stTRN,'LineWidth',2);
plot([-0.007:0.001:0.043-0.001],ratioProb1stS1,'LineWidth',2); 
plot([-0.007:0.001:0.043-0.001],ratioProb1stMOCK,'LineWidth',2); hold off
legend("VPM","TRN","S1","MOCK-VPM");
saveFigure(figure(1),['Z:\Jesus\LTP_Jesus_Emilio\PopFigures\TRN_VPM_S1figures/Ratio_ALL']);
figure;
plot([-0.007:0.001:0.043-0.001],ratioProb1stVPM,'LineWidth',2); hold on
plot([-0.007:0.001:0.043-0.001],ratioProb1stMOCK,'LineWidth',2);hold off
legend("VPM","MOCK-VPM"); title("Probability ratio After/Control");xlabel("time (s)")
ylabel("% Change")

saveFigure(figure(1),['Z:\Jesus\LTP_Jesus_Emilio\PopFigures\TRN_VPM_S1figures/Ratio_MockvsVPM'])
%%
CircuitModeAfter= NaN(300,3); CircuitModeAfter(1:size(VPMMode1stAfter),1)=VPMMode1stAfter;
CircuitModeAfter(1:size(TRNMode1stAfter),2)=TRNMode1stAfter; CircuitModeAfter(1:size(s1Mode1stAfter),3)=s1Mode1stAfter;
subplot(9,1,[1 2 3 4]); ax1=subplot(9,1,[1 2 3 4])
boxplot(CircuitModeAfter,'Orientation','horizontal');xlim([-0.010 0.043])
subplot(9,1,[6 7 8 9]); ax2=subplot(9,1,[6 7 8 9])
boxplot(CircuitMode,'Orientation','horizontal');xlim([-0.010 0.043])
linkaxes([ax1 ax2],'xy')

subplot(16,1,[5 6 7 8]); ax2=subplot(16,1,[5 6 7 8])
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stVPM(:),100)),'LineWidth',2);hold on; plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stTRN(:),100)),'LineWidth',2)
plot([-0.007:0.0005:0.043-0.0005],zscore(hist(ClustersL200Onset1stS1(:),100)),'LineWidth',2);hold off
xlim([-0.005 0.04])

subplot(16,1,[10 11 12 13 14 15 16]); ax3=subplot(16,1,[10 11 12 13 14 15 16])
secuantialactivationVPMTRNS1AFTERmap=cat(2,ClustersL200Onset1stVPM(:,4:45),ClustersL200Onset1stTRN(:,1:50),ClustersL200Onset1stS1(:,1:41));
NtotalCircuit=1:size(secuantialactivationVPMTRNS1AFTERmap,2);

imagesc(NtotalCircuit,[-0.007:0.001:0.043-0.001],zscore(hist(secuantialactivationVPMTRNS1AFTERmap,100)))
colormap(jet);caxis([-1 8]);view([-90 90]);hold on; plot([0;size(NtotalCircuit,2)],[0.0;0.0],'LineStyle','--','LineWidth',2,'Color','w')
hold off; ylim([-0.005 0.04])
linkaxes([ax1,ax2,ax3],'x');

