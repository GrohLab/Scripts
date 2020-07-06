% Choosing the working directory
dataDir = uigetdir('E:\Data\VPM\LTP','Choose a working directory');
if dataDir == 0
return
end
%dataDir = 'E:\Data\VPM\LTP\191016_Jesus_LTP_3710_1520_1500';
% dataDir = 'D:\LTP\190716_Jesus_Emilio LTP_3751_1520_1500';
figureDir = fullfile(dataDir,'Figures\');
load('D:\Ross\4.12.19\M162_Nt1CreChR2_all_channels.mat')
load('Z:\Ross\Experiments\smrx Files\4.12.19\M162_Nt1CreChR2_CFA_Sens_5ITI_CondSig.mat')
npoints5 = head68.npoints;
mechChan5 = chan68;
mObj5 = StepWaveform(chan68, fs);
load('Z:\Ross\Experiments\smrx Files\4.12.19\M162_Nt1CreChR2_CFA_Sens_10ITI_CondSig.mat')
mechChan10 = chan68;
npoints10 = head68.npoints;
mObj10 = StepWaveform(chan68, fs);
load('Z:\Ross\Experiments\smrx Files\4.12.19\M162_Nt1CreChR2_CFA_Sens_20ITI_CondSig.mat')
mechChan20 = chan68;
npoints20 = head68.npoints;
mObj20 = StepWaveform(chan68, fs);
Triggers5 = mObj5.subTriggers; offset5 = mObj5.NSamples;
Triggers10 = mObj10.subTriggers + offset5; offset10 = mObj10.NSamples;
Triggers20 = mObj20.subTriggers + offset5 + offset10;
TriggersAll = [Triggers5; Triggers10; Triggers20];
Triggers = TriggersAll;
npointsAll = npoints5 + npoints10 + npoints20;
Ns = npointsAll;
Nt = Ns/fs;
clc
load('D:\Ross\4.12.19\Adaptationanalysis.mat')
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
'UniformOutput',false);
Ncl = numel(goods);
mSubs = Triggers;

condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
MechTTL = mObj5.subs2idx(mSubs,Ns);
continuousSignals = {double(MechTTL)};
condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
[biGaps, gapSubs] = sort(diff(Conditions(chCond).Triggers(:,1)), 'descend');
gapFlag = abs(zscore(biGaps)) > 3;
mainGaps = sort(gapSubs(gapFlag), 'ascend');
Ng = numel(mainGaps);
Ncond = Ng+1;
condFlags = false(NTa, Ncond);
initSub = 1;
for cg = 1:Ng
condFlags(initSub:mainGaps(cg),cg) = true;
initSub = mainGaps(cg) + 1;
end
condFlags(initSub:NTa, Ncond) = true;
Na = sum(condFlags,1);
figure;imagesc(condFlags);




%% Counting spikes in given windows and computing the statistical significance
% Time logical flags
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
timeFlags = [sponActStackIdx;respActStackIdx];
% Time window
delta_t = diff(responseWindow);
% Statistical tests
[Results, Counts] = statTests(discStack,condFlags,timeFlags);
% Firing rate for all clusters, for all trials% Plotting
figs = scatterSignificance(Results, Counts, consCondNames, delta_t, sortedData(goods,1));
% meanfr = cellfun(@(x) mean(x,2)/delta_t,Counts,'UniformOutput',false);