%% light triggered psths

clear all
dataDir = 'Z:\Melina\Time Axis MatlabProblem- StartingPSTH.m';
dataDir = 'C:\Users\Rebecca Mease\Seafile\Time Axis MatlabProblem- StartingPSTH.m';


cd(dataDir)
binFiles = dir('*.bin'); %only works for some reason if we are in directory already
[~,expName,~] = fileparts(binFiles.name);  %a basic string for the experiment name RAM
expSubfix = fullfile(dataDir,expName);
%normally load data...
load board_dig_in_data
load ExportParameters fs
laserSignal = board_dig_in_data(3,:);%cheat for this one recording

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

% To Do: put this into a function that allows you to select a directory
% Name: loadExperiment

%% Getting the condition trigger
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

lstrength = [70,50,10,0,70,50,10,0]; %strength of laser
lsub=1; %which laser strength do you want to look at ?
nstim=100;  %How many stimuli were applied with this laserintensity?

%lstrength = [0,10,50,70]; %strength of laser
%lsub=5; %which laser strength do you want to look at ?
%nstim=100;  %How many stimuli were applied with this laserintensity?lstrength = [0,10,50,70]; %strength of laser

%lstrength = [70]; %strength of laser
%lsub=1; %which laser strength do you want to look at ?
%nstim=30;
subset = [1:nstim]+nstim*(lsub-1);
ltOn = ltOn(subset);
%ltOn=ltOn(1:100);  %YOU NEED TO CHANGE THIS FOR EVERY EXPERIMENT  (eg only
%considering the first 100 triggers)
%ltOn=ltOn(101:end);  %YOU NEED TO CHANGE THIS FOR EVERY EXPERIMENT

conditionString='30 trials'
%conditionString='laser zero'

clearvars lt
fprintf(1,'Got the condition triggers\n')

% To Do: put this into a function that extracts and saves trigger
% information
% Name: getTriggers


%%
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

sortedData(goods,:)  %display which units will be used
%% inferring neuron identify
% 1. channel information: user provides likely location based on channel
% map.  This need to be loaded from phy and is cluster_info.tsv
% A function for that takes channel groups and putative types and returns a
% cell array of sorted clusters and labels
%2. 

labels={'BC','POm'}  
regionAssignment={[207,407]... %cluster ID from channels in 'BC'  %#17; 50x_merged
   % [71,72,81,90,107,118,148,159,165,173,181,217,244,301,323,354,360,366,381,386,397]}  %cluster ID from channels in 'POm'
      [79,80,124,138,189,213,261,277,294,315,330,403]}        %VPM 
 
unit_labels=str2num(char(sortedData(goods,1)));

Groups={}
for i=1:numel(regionAssignment)
    [C,IA,IB]=intersect(regionAssignment{i},unit_labels)
    Groups{i}=goods(IB);
end

desiredUnits=Groups{1}  % this gives BC clusters
%desiredUnits=Groups{2}  % this gives POm clusters
%desiredUnits=goods
%%

% Logical spike trace for the first good cluster  % I don´t completely
% understand why E is doing this but I trust it RAM
spkLog = StepWaveform.subs2idx(round(sortedData{desiredUnits(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(desiredUnits(2:end),2),...
    'UniformOutput',false);
% Number of good clusters
Ncl = numel(goods);

% To Do: put this into a function that identifies and SAVES single units =
% 1
% and MUA= 2 and kicks out noise = 3  and take in argument for minimum
% firing rate
% Name: sortUnits
% this is called before trigger finding

%%


% this is somewhat redundant but oh well for now, set up for psth
lObj = StepWaveform(laserSignal,fs);
lSubs = lObj.subTriggers;
laser = lObj.subs2idx(lSubs,lObj.NSamples);
lObj.delete;
continuousSignals = {laser};
clearvars *Obj laser

%% User controlling variables to change details of PSTHs
% Time window to see the cluster activation in seconds

timeLapse = [0.02, 0.05];
%timeLapse = [0.040, 0.040];
%timeLapse = [0.040, 0.5];  % time before and time after, seconds (this variable will control output data that we save below)
%timeLapse = [0.01, 0.7];

% Bin size for PSTHs
binSz = 0.0005; %in seconds
%binSz = 0.0005; % 0,5ms
consideredConditions = 1;
Nccond = length(consideredConditions);
% Adding all the triggers from the piezo and the laser in one array
allWhiskersPlusLaserControl = ltOn; %ugh
% Logical and numerical stack for computations
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

[PSTH, trig, sweeps] = getPSTH(dst,timeLapse,false(size(ltOn)),binSz,fs);

fig = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
    [{'Laser'};sortedData(desiredUnits,1)],conditionString);
 axes(fig.Children(2))
 x=get(gca,'YAxis')
%x(1).Limits=[0 .01]   %% change the axis
x(1).Limits=[0 .2]   %% change the axis

%axes(fig.Children(3))  %decoment this block to see the colorplots in blue!
%colormap parula
%shading interp

configureFigureToPDF(fig);
print(fig,fullfile(dataDir,sprintf('%s %s.pdf',...
    expName, conditionString)),...
    '-dpdf','-fillpage')

% configureFigureToPDF(fig);
% print(fig,fullfile(dataDir,sprintf('%s %s.eps',...
%     expName, 'laser')),...
%     '-deps')

figureFileName=fullfile(dataDir,sprintf('%s %s.fig',expName, conditionString))
savefig(fig,figureFileName)   %saves population psths

%think about automating 

%
ClusterIds=sortedData(desiredUnits,1);
tx = linspace(-timeLapse(1),timeLapse(2),Nt);
psthTX = linspace(-timeLapse(1),timeLapse(2),size(PSTH,2));

Pop={};
for i=2:size(dst,1) %for every unit
    Sp={};
     for j=1:size(dst,3) %for every trial
         Sp{j}=tx(find(squeeze(dst(i,:,j))));
     end
     Pop{i-1}=cell2mat(Sp);
end


%%

figure
nsp=floor(sqrt(numel(Pop)));

 %overall scaling for y axis
%look here tyscaling = 1.05 o change parameters if too many neurons
tiled=true
tall=1% column
wide=2 %row
count=0
for i=1:numel(Pop)
    if numel(Pop{i})>10
        count=count+1;
        if tiled
            subplot(tall,wide,count)
        else
            figure
        end
        h=histogram(Pop{i}*1000,'binwidth',binSz*1000,'Normalization','count')
        %xtick_ms=10.^get(gca,'XTick');
        %set(gca,'XTickLabel',xtick_ms)
        title([ClusterIds{i} '_ClusterId'], 'Interpreter', 'none')
        xlabel('ms')
        ylabel('event counts')
        hold on
        plot(psthTX*1000, PSTH(i,:))
        axis tight
        ylim([0 max(max(PSTH))*yscaling])
    else
    end
end

%% now manual change if necessary
%  click on figure, then run block of code
lims=[100 500] %manually insert limits here
count=0
numplots=2

tall=1% column
wide=2 %row

for i=1:numplots
    subplot(tall,wide,i)
    ylim([0 lims(i)])
end


%%
%if you wanted to zoom
for i=1:numplots
    subplot(tall,wide,i)
    xlim([-1 10])
end


%%
print -dmeta
print(gcf,'allPSTHs','-dmeta')

%% Save an image: click on this section while the figure to save is open then it gets saved as PDF
fig = gcf;
name = get(get(gca,'Title'),'String');

configureFigureToPDF(fig);
print(fig,fullfile(dataDir,sprintf('%s %s.pdf',expName,[conditionString name])),'-dpdf','-fillpage')

%%
%allSelectionIdx = true(numel(goods),1);
%PopRelativeSpikeTimes=getRasterFromStack(dst,false(size(ltOn)),allSelectionIdx, timeLapse, fs, true);  %



%% 

%close all
figure
Spikes=sortedData(desiredUnits,2);

binsize=.050 %in seconds CHANGE THIS AND LOOK IN PCA
%binsize=0.0001
spikes =Spikes{:};
H=histogram(spikes,'BinWidth',binsize)   

BinEdges=H.BinEdges;
Rs=[];

%change these if need be=========================
tiled=true
tall=4 %rows
wide=5 %columns
maxSubPanels=20; 
%=========================

count=0
for i=1:numel(Spikes)
    spikes=Spikes{i};
    if count==0
        figure
    end
    if tiled
        count=count+1;
        subplot(tall,wide,count)
    end
    
    %make new fig if too many CHANGE THISd
    if count>maxSubPanels
        count=0;
    end
    
    H=histogram(spikes,BinEdges)
    title([ 'Id=' ClusterIds{i}  ' #=' num2str(i)])
    xlabel('time in sec')
    ylabel('spike rate in Hz') %spike counts normalized and bin size is in seconds
    R=H.Values;
    Rs=[Rs R'];
end


%TO DO-- use locfit from Chronux package to get smoother versions of r(t) 

%% correlation between rates
[rcoeff p]=corrcoef(Rs);
figure
rcoeff(find(rcoeff==1))=nan;
subplot(1,2,1)
pcolor(rcoeff), shading flat
colorbar
subplot(1,2,2)
hist(rcoeff(:))


[i,j] = find(p<0.05);  % find significant correlations
[i,j]                  % display their (row,col) indices

[i,j] = find(rcoeff>0.5);  % find significant correlations

[i,j,rcoeff(ind2sub(size(rcoeff),find(rcoeff>0.5)))]             % display their (row,col) indices


%%
   [COEFF, SCORE, LATENT] = pca(Rs)
   figure
   plot(1:numel(LATENT),cumsum(LATENT)/sum(LATENT))
   
   %%
   figure
   plot(SCORE(:,1:2))
    
    
    %%
    
    figure
    plot(SCORE(:,1), SCORE(:,2),'.')
    figure
    plot3(SCORE(:,1), SCORE(:,2), SCORE(:,3),'.')
    %%
    figure
    scatter3(SCORE(:,1), SCORE(:,2), SCORE(:,3),15,(SCORE(:,4)),'filled')
