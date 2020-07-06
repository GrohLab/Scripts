% 28.11.19 anesthetized rbp4 opto whisker analysis by using the Data Explorer. 
clearvars
%% Load the data
%dataDir = 'F:\Kilosorting\MC-1-2019\191018_anesthetized\Combined recording_kilosorting again';
dataDir = 'F:\Kilosorting\MC-1-2019\191018_anesthetized\Concatination single files\Laser';
cd(dataDir)
binFiles = dir('*.bin'); %only works for some reason if we are in directory already
[~,expName,~] = fileparts(binFiles.name);  %a basic string for the experiment name RAM

% Loading the sampling frequency, the sorted clusters, and the conditions
% and triggers.
expSubfix = fullfile(dataDir,expName);
try
    load([expSubfix,'_sampling_frequency.mat'],'fs')
catch
    fprintf(1,'Seems like the kilosort-phy pipeline hasn''t been touched!\n')
    return
end

%get spikes
try
    load([expSubfix,'_all_channels.mat'],'sortedData')
catch
    try
        %importPhyFiles(dataDir,expName,false,true);  %we will try to load channel info as well
        importPhyFiles(dataDir)
    catch
        fprintf(1,'Error importing the phy files into Matlab format\n')
        return
    end
    load([expSubfix,'_all_channels.mat'],'sortedData')  %this needs to be changed to load in channel info as well
end

%%  get trigger information
%load analysis files if it exists

if exist([expSubfix, '_CondSig.mat'],'file')   %if CondSig exists, then directly call delay protocol
    getDelayProtocol(dataDir);
else   %if not, then make it
    getConditionSignalsBF(fopen([expSubfix,'.smrx']))
    getDelayProtocol(dataDir);
end
load([expSubfix, '_CondSig.mat']) %loading the triggers

laserSignal=chan35;
whiskerSignal=chan34;

%% Getting the condition triggers
% These variables detect the rising and falling edges of the digital input
% pulse (laser).

stObj = StepWaveform(laserSignal,fs,'On/Off','Laser');
%clearvars laserSignal %keep this for later
% The lt variable contains a logical array size 2(on/off) x
% total_number_of_samples indicating the beginning and end of a pulse with
% a true value.
lt = stObj.Triggers;
% This contains the sample number at which the pulse rose. If you want to
% change for pulse fall, write as follows: ltOn = find(lt(:,2));
ltOn = find(lt(:,1));

 %IF YOU WANT TO LOOK AT ALIGNMENT OF DIG IN AND TRIGGERS
figure
plot(laserSignal)
hold on
plot(ltOn,laserSignal(ltOn),'o')

% These variables detect the rising and falling edges of the digital input
% pulse (whisker).
stObj = StepWaveform(whiskerSignal,fs,'On/Off','Whisker');
%clearvars laserSignal %keep this for later
% The lt variable contains a logical array size 2(on/off) x
% total_number_of_samples indicating the beginning and end of a pulse with
% a true value.
wt = stObj.Triggers;
% This contains the sample number at which the pulse rose. If you want to
% change for pulse fall, write as follows: ltOn = find(lt(:,2));
wOn = find(wt(:,1));

 %IF YOU WANT TO LOOK AT ALIGNMENT OF DIG IN AND TRIGGERS
hold on
plot(whiskerSignal)
hold on
plot(wOn,whiskerSignal(wOn),'o')

%ltOn=wOn   %GET RID OF THIS TO LOOK AT LASER   
lstrength = [100]; %strength of laser
lsub=1; %which laser strength do you want to look at ?
nstim=50;
subset = [1:nstim]+nstim*(lsub-1);
ltOn = ltOn(subset);
%ltOn=ltOn(1:100);  %YOU NEED TO CHANGE THIS FOR EVERY EXPERIMENT  (eg only
%considering the first 100 triggers)
%ltOn=ltOn(101:end);  %YOU NEED TO CHANGE THIS FOR EVERY EXPERIMENT

conditionString='50 trials laser'
%conditionString='laser zero'

clearvars lt
fprintf(1,'Got the condition triggers\n')

% grabbing good spikes
Ns = numel(laserSignal)
% Total duration of the recording
Nt =Ns/fs;  %seconds

% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2)); %count spikes per unit
clusterSpikeRate = totSpkCount/Nt;  %divide by seconds to get rate
MinimumSpikeRate=0.01;  % in Hz change this to vary minimum spike rate allowed
silentUnits = clusterSpikeRate < MinimumSpikeRate;  %find low rate units
bads = union(bads,find(silentUnits)); %add low rate units to bad category (merging real neural signal with noise!)
goods = setdiff(1:size(sortedData,1),bads); %define the good units as not bad
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
%
% Logical spike trace for the first good cluster  % I don´t completely
% understand why E is doing this but I trust it RAM
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
    'UniformOutput',false);
% Number of good clusters
Ncl = numel(goods);


% this is somewhat redundant but oh well for now
lObj = StepWaveform(laserSignal,fs);
lSubs = lObj.subTriggers;
laser = lObj.subs2idx(lSubs,lObj.NSamples);
lObj.delete;
continuousSignals = {laser};
clearvars *Obj laser

%% User controlling variables to change details of PSTHs
% Time window to see the cluster activation in seconds
timeLapse = [0.050, 0.250];  % time before and time after, seconds (this variable will control output data that we save below)
% Bin size for PSTHs
binSz = 0.0005; %in seconds
consideredConditions = 1;
Nccond = length(consideredConditions);
% Adding all the triggers from the piezo and the laser in one array
allWhiskersPlusLaserControl = ltOn; %ugh
%% Logical and numerical stack for computations
% dst - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.
[dst, cst] = getStacks(spkLog, ltOn,'on',timeLapse,fs,fs,spkSubs,continuousSignals);
% Number of clusters + the piezo as the first event + the laser as the last
% event, number of time samples in between the time window, and number of
% total triggers.
[Ne, Nt, NTa] = size(dst);   %got confused because it meant something else earlier RAM
% Computing the time axis for the stack
tx = (0:Nt)/fs - timeLapse(1);

binSize=.050  %in ms
[PSTH, trig, sweeps] = getPSTH(dst,timeLapse,false(size(ltOn)),binSz,fs);

fig = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
    [{'Laser'};sortedData(goods,1)],conditionString);

configureFigureToPDF(fig);
print(fig,fullfile(dataDir,sprintf('%s %s.pdf',...
    expName, conditionString)),...
    '-dpdf','-fillpage')

% configureFigureToPDF(fig);
% print(fig,fullfile(dataDir,sprintf('%s %s.eps',...
%     expName, 'laser')),...
%     '-deps')
savefig(fig, fullfile(dataDir,sprintf('%s %s.fig',expName, conditionString)))   %saves population psths



%%
ClusterIds=sortedData(goods,1);
tx = linspace(-timeLapse(1),timeLapse(2),Nt);

Pop={};
for i=2:size(dst,1) %for every unit
    Sp={};
     for j=1:size(dst,3) %for every trial
         Sp{j}=tx(find(squeeze(dst(i,:,j))));
     end
     Pop{i-1}=cell2mat(Sp);
 end
 
for i=1:numel(Pop)
     if numel(Pop{i})>10
        figure
        histogram(Pop{i}*1000,200)
        %xtick_ms=10.^get(gca,'XTick');
        % set(gca,'XTickLabel',xtick_ms)
        title([ClusterIds{i} '_ClusterId'], 'Interpreter', 'none')
        xlabel('ms')
        ylabel('event counts')
        xlim([-100 150])
    else
     end
end



%%  look at spontaneous activity
load SpikingParameters
% just look at the spike times for one/all cluster
%Visualize spikes over time --> raster plot for a given chunk of time

close all
minTime=0 %in seconds
maxTime=4*60  %in seconds  %check that this is within the range of the recording

%spikes=Spikes{1}; %grab spike times from first cluster
figure
level=0

for i=1:numel(Spikes)  %grab spike times from all "good" cluster
    level=level+1;
    spikes=Spikes{i};
    spikes=spikes(spikes>minTime & spikes<maxTime);
    y=ones(size(spikes))*level;
    plot(spikes,y,'LineStyle','none','Marker','.','Markersize',10)
  
    ylabel('clusters')
    hold on
end
 xlabel('time in sec')