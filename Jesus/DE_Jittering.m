% 30.08.19 Jittering analysis by using the Data Explorer. 
clearvars
%% Load the data
% Choosing the working directory
dataDir = uigetdir('E:\Data\VPM\LTP','Choose a working directory');
if dataDir == 0
    return
end
%dataDir = 'E:\Data\VPM\LTP\191016_Jesus_LTP_3710_1520_1500';
% dataDir = 'D:\LTP\190716_Jesus_Emilio LTP_3751_1520_1500';
figureDir = fullfile(dataDir,'Figures\');
% Loading the necessary files
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables\n')
    return
end
%{
% dataDir = 'E:\Data\VPM\LTP\190703_LTP_3720_1520_1520\LTP2';
% binFiles = dir(fullfile(dataDir,'*.bin'));
% [~,expName,~] = fileparts( binFiles.name);
% figureDir = fullfile(dataDir, 'Figures');
% % Loading the sampling frequency, the sorted clusters, and the conditions
% % and triggers.
% expSubfix = fullfile(dataDir,expName);
% try
%     load([expSubfix,'_sampling_frequency.mat'],'fs')
% catch
%     fprintf(1,'Seems like the kilosort-phy pipeline hasn''t been touched!\n')
%     return
% end
% try
%     load([expSubfix,'analysis.mat'],'Conditions','Triggers')
% catch
%     if exist([expSubfix, '_CondSig.mat'],'file')
%         getDelayProtocol(dataDir);
%     else
%         try
%             getConditionSignalsBF(fopen([expSubfix,'.smrx']))
%             getDelayProtocol(dataDir);
%         catch
%             fprintf(1,'Confusing naming. Cannot continue\n')
%             return
%         end
%     end
%     load([expSubfix,'analysis.mat'],'Conditions','Triggers')
% end
% try
%     load([expSubfix,'_all_channels.mat'],'sortedData')
% catch
%     try
%         importPhyFiles(dataDir);
%     catch
%         fprintf(1,'Error importing the phy files into Matlab format\n')
%         return
%     end
%     load([expSubfix,'_all_channels.mat'],'sortedData')
% end
%}
%% Constructing the helper 'global' variables
% Number of total samples
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;
% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
    'UniformOutput',false);
% Number of good clusters 
Ncl = numel(goods);
% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker'};
cxStim = {'laser','light'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
while ctn <= numTrigNames 
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        whisker = Triggers.(trigNames{ctn});
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        laser = Triggers.(trigNames{ctn});
    end
    ctn = ctn + 1;
end
mObj = StepWaveform(whisker,fs);
mSubs = mObj.subTriggers;
piezo = mObj.subs2idx(mSubs,Ns);
lObj = StepWaveform(laser,fs);
lSubs = lObj.subTriggers;
laser = lObj.subs2idx(lSubs,Ns);
mObj.delete;lObj.delete;
continuousSignals = {piezo;laser};
clearvars *Obj piezo laser
%% User controlling variables
%{
% Time window to see the cluster activation in seconds
timeLapse = [0.05, 0.15];
% Bin size for PSTHs
binSz = 0.0005;
%}
% Time lapse, bin size, and spontaneous and response windows
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'0.1, 0.1', '0.002, 0.05', '0.001'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n')
    return
else
    timeLapse = str2num(answ{1}); %#ok<*ST2NM>
    if numel(timeLapse) ~= 2
        timeLapse = str2num(inputdlg('Please provide the time window [s]:',...
            'Time window',[1, 30], '0.1, 0.1'));
        if isnan(timeLapse) || isempty(timeLapse)
            fprintf(1,'Cancelling...')
            return
        end
    end
    responseWindow = str2num(answ{2});
    binSz = str2double(answ(3));
end
fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
spontaneousWindow = -flip(responseWindow);
%% Condition triggered stacks
condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
condGuess = contains(condNames, 'whiskerall', 'IgnoreCase', true);
% Choose the conditions to create the stack upon
[chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
    'PromptString',...
    'Choose the condition which has all whisker triggers: (one condition)',...
    'InitialValue', find(condGuess), 'ListSize', [350, numel(condNames)*16]);
if ~iOk
    fprintf(1,'Cancelling...\n')
    return
end

% Select the onset or the offset of a trigger
fprintf(1,'Condition ''%s''\n', Conditions(chCond).name)
onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
    'on','off','Cancel','on');
if strcmpi(onOffStr,'Cancel')
    fprintf(1,'Cancelling...\n')
    return
end

% Constructing the stack out of the user's choice
% discStack - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.
[discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);

% [dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
%     'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
%     continuousSignals);

% Number of clusters + the piezo as the first event + the laser as the last
% event, number of time samples in between the time window, and number of
% total triggers.
[Ne, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs - timeLapse(1);
%% Considered conditions selection
% Choose the conditions to look at
auxSubs = setdiff(1:numel(condNames), chCond);
ccondNames = condNames(auxSubs);
[cchCond, iOk] = listdlg('ListString',ccondNames,'SelectionMode','multiple',...
    'PromptString',...
    'Choose the condition(s) to look at (including whiskers):',...
    'ListSize', [350, numel(condNames)*16]);
if ~iOk
    fprintf(1,'Cancelling...\n')
    return
end

% Select the onset or the offset of a trigger
fprintf(1,'Condition(s):\n')
fprintf('- ''%s''\n', Conditions(auxSubs(cchCond)).name)


% Subscript to indicate the conditions with all whisker stimulations,
% whisker control, laser control, and the combination whisker and laser.
allWhiskerStimulus = chCond;
% allLaserStimulus = 2;
% whiskerControl = 9;
% laserControl = 8;
consideredConditions = auxSubs(cchCond);

Nccond = length(consideredConditions);
% Time windows to evaluate if a unit is responsive or not.
% spontaneousWindow = [-0.05, -0.002];
% responseWindow = [0.002, 0.05];
% Adding all the triggers from the piezo and the laser in one array
% allWhiskersPlusLaserControl = ...
%     union(Conditions(allWhiskerStimulus).Triggers,...
%     Conditions(laserControl).Triggers,'rows');
%% Boolean flags 
delayFlags = false(NTa,Nccond);
counter2 = 1;
for ccond = consideredConditions
    delayFlags(:,counter2) = ismember(Conditions(chCond).Triggers(:,1),...
        Conditions(ccond).Triggers(:,1));
    counter2 = counter2 + 1;
end
Na = sum(delayFlags,1);
%% Computing which units/clusters/putative neurons respond to the stimulus
% Logical indices for fetching the stack values
sponActStackIdx = tx >= spontaneousWindow(1) & tx <= spontaneousWindow(2);
respActStackIdx = tx >= responseWindow(1) & tx <= responseWindow(2);
% The spontaneous activity of all the clusters, which are allocated from
% the second until one before the last row, during the defined spontaneous
% time window, and the whisker control condition. 

timeFlags = [sponActStackIdx;respActStackIdx];
% Time window
delta_t = diff(responseWindow);
% Statistical tests
[Results, Counts] = statTests(discStack,delayFlags,timeFlags);
% Firing rate for all clusters, for all trials
meanfr = cellfun(@(x) mean(x,2)/delta_t,Counts,'UniformOutput',false);

indCondSubs = cumsum(Nccond:-1:1);
consCondNames = condNames(consideredConditions);


%{
sponTimeMarginal = sum(...
    discStack(2:Ne-1,sponActStackIdx,delayFlags(:,Nccond)),2);
sponTimeMarginal = squeeze(sponTimeMarginal);
sponActPerTrial = sum(sponTimeMarginal,2)/Na(Nccond);
% Similarly for the responsive user-defined time window
respTimeMarginal = sum(...
    discStack(2:Ne-1,respActStackIdx,delayFlags(:,Nccond)),2);
respTimeMarginal = squeeze(respTimeMarginal);
respActPerTrial = sum(respTimeMarginal,2)/Na(Nccond);
activationIndex = -log(sponActPerTrial./respActPerTrial);
%}
whiskerResponsiveUnitsIdx = activationIndex > 1;
display(find(whiskerResponsiveUnitsIdx))
%% Getting the relative spike times for the whisker responsive units (wru)
% For each condition, the first spike of each wru will be used to compute
% the standard deviation of it.
cellLogicalIndexing = @(x,idx) x(idx);
isWithinResponsiveWindow =...
    @(x) x > responseWindow(1) & x < responseWindow(2);

Nwru = sum(whiskerResponsiveUnitsIdx);
unitSelectionIdx = [whiskerResponsiveUnitsIdx(2:end);false];
firstSpike = zeros(Nwru,Nccond);

for ccond = 1:size(delayFlags,2)
    relativeSpikeTimes = getRasterFromStack(discStack,~delayFlags(:,ccond),...
        unitSelectionIdx, timeLapse, fs, true, true);
    relativeSpikeTimes(~whiskerResponsiveUnitsIdx(1),:) = [];
    respIdx = cellfun(isWithinResponsiveWindow, relativeSpikeTimes,...
        'UniformOutput',false);
    spikeTimesINRespWin = cellfun(cellLogicalIndexing,...
        relativeSpikeTimes, respIdx, 'UniformOutput',false);
    for ccl = 1:Nwru
        frstSpikeFlag = ~cellfun(@isempty,spikeTimesINRespWin(ccl,:));
        firstSpike(ccl,ccond) = std(...
            cell2mat(spikeTimesINRespWin(ccl,frstSpikeFlag)));    
    end
end
% This line looks pretty for saving the relative spike times.
%% Plotting the population activity
% On the fly section: goodsIdx is the negated version of badsIdx
goodsIdx = ~badsIdx';
for ccond = 1:Nccond
    [PSTH, trig, sweeps] = getPSTH(... dst([true;whiskerResponsiveUnitsIdx;true],:,:),timeLapse,...
        discStack,timeLapse,...
        ~delayFlags(:,ccond),binSz,fs);
    fig = plotClusterReactivity(PSTH,trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};... sortedData(goods(whiskerResponsiveUnitsIdx),1);{'Laser'}],...
        sortedData(goods,1);{'Laser'}],...
        strrep(expName,'_','\_'));
    configureFigureToPDF(fig);
    print(fig,fullfile(dataDir,sprintf('%s %s.pdf',...
        expName, Conditions(consideredConditions(ccond)).name)),...
        '-dpdf','-fillpage')
    print(fig,fullfile(dataDir,sprintf('%s %s.emf',...
        expName, Conditions(consideredConditions(ccond)).name)),...
        '-dmeta')
end

