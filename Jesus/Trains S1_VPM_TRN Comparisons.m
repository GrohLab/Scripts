%% primer stimulo tren VPM TRN S1 Control

%ClustersControlOnsetVPMbad=ClustersControlOnsetVPM>=0 & ClustersControlOnsetVPM<=0.002
%ClustersControlOnsetVPMbad1=ClustersControlOnsetVPM>=0.101 & ClustersControlOnsetVPM>=0.104;
%ClustersControlOnsetVPMbad=(ClustersControlOnsetVPM>=0 & ClustersControlOnsetVPM<=2) | (ClustersControlOnsetVPM>=0.100&ClustersControlOnsetVPM<=0.104);
%ClustersControlOnsetVPM1=ClustersControlOnsetVPM(ClusetersControlONsetVPMbad~=1);
ClustersControlOnsetVPMms=(ClustersControlOnsetVPM*1000)+3;
ClustersControlOnsetTRNms=(ClustersControlOnsetTRN*1000)+3;
ClustersControlOnsetS1ms=(ClustersControlOnsetS1*1000)+3;
TimePsth=[-197:1:1203-1]
Binsize=0.001

ClustersControlOnsetVPMTren=ClustersControlOnsetVPMms(:,54:end);
NclustersVPMTren=size(ClustersControlOnsetVPMms(:,54:end),2);
NclustersTRNtren=size(ClustersControlOnsetTRNms,2);
NclustersS1tren=size(ClustersControlOnsetS1ms,2);

NormControlVPM=hist(ClustersControlOnsetVPMTren(:),1400)./(150*(NclustersVPMTren*Binsize));
figure(1);
plot(TimePsth,NormControlVPM)
hold on
NormControlTRN=hist(ClustersControlOnsetTRNms(:),1400)./(150*NclustersTRNtren*Binsize);
plot(TimePsth,NormControlTRN)
NormControlS1=hist(ClustersControlOnsetS1ms(:),1400)./(150*NclustersS1tren*Binsize)
plot(TimePsth,NormControlS1);
legend("VPM","TRN","S1Bf");
xlabel('time (ms)'); ylabel('Firing Rate (Hz)');title('Whisker system Psth Control');xlim([-150 1150]);

saveFigure(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/HzCircuitPsthControl']);

TimePsth2=[-200:0.5:1200-0.5];
%probabilidad VPM
HistVPM=hist(ClustersControlOnsetVPMTren(:),1400); 
ProbVPM=HistVPM./(sum(HistVPM));
%probabilidad TRN
HistTRN=hist(ClustersControlOnsetTRNms(:),1400);
ProbTRN=HistTRN./(sum(HistTRN));
%probability S1
HistS1=hist(ClustersControlOnsetS1ms(:),1400);
ProbS1=HistS1./(sum(HistS1)); figure(2);
plot(TimePsth,ProbVPM); hold on;
plot(TimePsth,ProbTRN);
plot(TimePsth,ProbS1); hold off;
legend("VPM","TRN","S1Bf");
xlabel('time (ms)'); ylabel('Probability');title('Probability Whisker system Psth Control'),;xlim([-150 1150]);
%save figure
%foldername='Z:\Jesus\LTP_Jesus_Emilio\PopFigures\TRN_VPM_S1figures'
%saveas(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/ProbCircuitPsth.emf']);
%saveas(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/ProbCircuitPsth.pdf']);
%saveas(figure(1),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/ProbCircuitPsth.fig']);

saveFigure(figure(2),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/ProbCircuitPsthControl']);
%% afterinduction
ClustersL200OnsetVPMms=(ClustersL200OnsetVPM*1000)+3;
ClustersL200OnsetTRNms=(ClustersL200OnsetTRN*1000)+3;
ClustersL200OnsetS1ms=(ClustersL200OnsetS1*1000)+3;
TimePsth=[-197:1:1203-1]
Binsize=0.001

ClustersL200OnsetVPMTren=ClustersL200OnsetVPMms(:,54:end);
NclustersVPMTrenafter=size(ClustersL200OnsetVPMms(:,54:end),2);
NclustersTRNtrenafter=size(ClustersL200OnsetTRNms,2);
NclustersS1trenafter=size(ClustersL200OnsetS1ms,2);
figure;
NormControlVPMafter=hist(ClustersL200OnsetVPMTren(:),1400)./(300*NclustersVPMTrenafter*Binsize);
plot(TimePsth,NormControlVPMafter)
hold on
NormControlTRNafter=hist(ClustersL200OnsetTRNms(:),1400)./(300*NclustersTRNtrenafter*Binsize);
plot(TimePsth,NormControlTRNafter)
NormControlS1after=hist(ClustersL200OnsetS1ms(:),1400)./(300*NclustersS1trenafter*Binsize)
plot(TimePsth,NormControlS1after); hold off
legend("VPM","TRN","S1Bf");
xlabel('time (ms)'); ylabel('Firing Rate (Hz)');title('Whisker system Psth After-Induction');;xlim([-150 1150]);

saveFigure(figure,['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/HzCircuitPsth-After']);

%probabilidad VPM
HistVPMAfter=hist(ClustersL200OnsetVPMTren(:),1400); 
ProbVPMAfter=HistVPMAfter./(sum(HistVPMAfter));
%probabilidad TRN
HistTRNAfter=hist(ClustersL200OnsetTRNms(:),1400);
ProbTRNAfter=HistTRNAfter./(sum(HistTRNAfter));
%probability S1
HistS1After=hist(ClustersL200OnsetS1ms(:),1400);
ProbS1After=HistS1After./(sum(HistS1After)); figure(2);
plot(TimePsth,ProbVPMAfter); hold on;
plot(TimePsth,ProbTRNAfter);
plot(TimePsth,ProbS1After); hold off;
legend("VPM","TRN","S1Bf");
xlabel('time (ms)'); ylabel('Probability');title('probability Whisker system Psth After-Inductiom');;xlim([-150 1150]);
%save figure
%foldername='Z:\Jesus\LTP_Jesus_Emilio\PopFigures\TRN_VPM_S1figures'

saveFigure(figure(2),['Z:/Jesus\LTP_Jesus_Emilio/PopFigures/TRN_VPM_S1figures/Prob_CircuitPsth-After']);