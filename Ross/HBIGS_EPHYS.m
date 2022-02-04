% HBIGs 2022 in vivo electrophysiology analysis tutorial

% Welcome to the HBIGS tutorial. In this tutorial will be be covering the
% basics of the kinds of analysis we do in the Groh Lab, specifically
% multi-unit extracellular electrophysiology.

% Whilst the script should run in full, we recommend (especially for those
% unfamiliar with MATLAB or programming) to run the script section by
% section to really get a feel for how and why each step is necessary for answering the questions that we're interested in.
% You can do this run each section by either clicking th 'Run Section' icon in the
% Editor section of the taskbar, or by pressing Strg + Enter.
% To run a section and then automatically advance to the next section in
% the script there is also a 'Run and Advance' button, or alternatively
% Strg + Shift + Enter will do this.



%% Firstly, we need a data directory to work in:

dataDir = 'Z:\HBIGS';

% We can see what files are in the current directory with the ls functino:

ls(dataDir)

%% Concatenating and converting smrx files to bin file.

% The .smrx files contain the data that was collected during the experiment.
% Therefore we need to access both the voltage traces from the probes, and
% therefore the traces of any non-neural paramters we may be interested in
% (e.g. TTL traces, respiration rate, etc.)


% To take votlage traces from multiple electrodes and and assign the
% spikes from these into respective neurons, we need to put these traces
% through a 'Spike-sorter'; an algorithm that clusters similar spikes as
% belonging to these putative neurons (referred to as units). But first,
% we need to concatenate the separate smrx files together and convertthem into a single binary file to be read by the
% spike-sorting software.


msmrx2bin('Z:\HBIGS','HBIGS_Tutorial')

% The chronological order of the recording are denoted by 'P', where 'P1_'
% denotes the first recording of the experiment.



%% Now that we have the bin file, we can run our spike-sorter.


cd C:\Users\NeuroNetz\Documents\GitHub\Kilosort2_5

kilosort

% probe file = 'Z:\Ross\ProbeFiles\' - Corrected_H3_ChanMap.mat
fs = 3.003003003003003e+04;




%% Once sorted and curated we can import our units and their spikes into MATLAB.

cd C:\Users\NeuroNetz\Documents\GitHub\Scripts\Ross

% Experiment Name
expName = 'I_hope_this_works_please!';

% Spike-sorted units and their spike times
sortedData = importPhyFiles('Z:\HBIGS\Pre-sorted', expName);

%Unit/Cluster (same thing) information in a table.
clInfo = getClusterInfo([dataDir, '\cluster_info.tsv']);


%% Now we need to organise the TTL data to know when we are stimulating the brain.

% Calling a function to give us the Conditions and the Trigger times.

% We will need the MechTTL, MechStim, and Laser channels
[Conditions, Triggers] = getTriggers(expName, dataDir, fs);



%% If any problems occur, we can also jut load the spikes and triggers that have been pre-sorted and orgainised earlier

load('Z:\PainData\Corrected_Channel_Map\L6\Cortex\20.8.21\KS2\m6_newChanMap_all_channels.mat');
load('Z:\PainData\Corrected_Channel_Map\L6\Cortex\20.8.21\KS2\m6_analysis.mat');
clInfo = getClusterInfo('Z:\PainData\Corrected_Channel_Map\L6\Cortex\20.8.21\KS2\cluster_info.tsv');
expName = 'Replace_this_text_with_your_experiment_name';
%% Creating a directory to put our figures in to.


% Creating the figure directory
figureDir = fullfile(dataDir,'Figures\');
if ~mkdir(figureDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end


%% Constructing the helper 'global' variables

% Triggers.MechStim = Triggers.MechStim * -1;

% Number of total samples
Ns = min(structfun(@numel,Triggers));

% Total duration of the recording in seconds
Nt = Ns/fs;

% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = badsIdx | silentUnits;

% Indicating which units will be analysed in the clInfo table
% (ActiveUnits!)
if ~any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
    clInfo = addvars(clInfo,~badsIdx,'After','id',...
        'NewVariableNames','ActiveUnit');
end

% Get the IDs of the axctive units
gclID = sortedData(goods,1);
badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));

% Logical spike trace for the first good cluster
spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);


% Subscript column vectors for the rest good clusters
spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
    'UniformOutput',false);

% Number of good clusters
Ncl = numel(goods);

% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker','mech','audio'};
cxStim = {'laser','light'};
lfpRec = {'lfp','s1','cortex','s1lfp'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
continuousSignals = cell(numTrigNames,1);
continuousNameSub = zeros(size(trigNames));
while ctn <= numTrigNames
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},lfpRec,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    ctn = ctn + 1;
end
continuousSignals(continuousNameSub == 0) = [];
continuousNameSub(continuousNameSub == 0) = [];
trigNames = trigNames(continuousNameSub);




%% Now we have identified the units we want to analyse, and the different manipulations we've induced in thenm (Conditions). 

% With multi-unit extracellular recordngs, the information we reliably have
% access to and can answer question about, are unit firing (spike) rates
% and the inter-spike intervals. 

% Lets calculate the interspike intervals for later

isiFile = fullfile(dataDir,[expName,'_ISIvars.mat']);
if ~exist(isiFile,'file')
    spkSubs2 = cellfun(@(x) round(x.*fs), sortedData(goods,2),...
        'UniformOutput', false);
    ISIVals = cellfun(@(x) [x(1)/fs; diff(x)/fs], spkSubs2,...
        'UniformOutput', 0);
    NnzvPcl = cellfun(@numel,ISIVals);
    Nnzv = sum(NnzvPcl);
    rows = cell2mat(arrayfun(@(x,y) repmat(x,y,1), (1:Ncl)', NnzvPcl,...
        'UniformOutput', 0));
    cols = cell2mat(spkSubs2);
    vals = cell2mat(ISIVals);
    ISIspar = sparse(rows, cols, vals);
else
    load(isiFile,'ISIspar')
end



%% We can now decide how to porceed with respect to the time around each
% manipulation that we want to consider. 



% Time lapse, bin size, and spontaneous and response windows

    % Viewing window denotes the timespan around the TTL trigger shown for the
    % Peri-Stimulus-Time-Histograms (PSTHs - something we'll get on to later)
    
    % Response window is the timespan after the TTL trigger to be
    % considered for the statistical analysis.
    
    % Bin size is also a PSTH measure that determines how many spikes are bunched into a given bin. 
    
    % Feel free to play around with these parameters once after first
    % running with the defaults.
    
    
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'-2, 6', '0.5, 2', '0.01'};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
if isempty(answ)
    fprintf(1,'Cancelling...\n');
    return
else
    timeLapse = str2num(answ{1}); %#ok<*ST2NM>
    if numel(timeLapse) ~= 2
        timeLapse = str2num(inputdlg('Please provide the time window [s]:',...
            'Time window',[1, 30], '-0.1, 0.1'));
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


% In most scenarios it makes sense to take the same time-window before each
% TTL trigger to compare activity statistically

sponAns = questdlg('Mirror the spontaneous window?','Spontaneous window',...
    'Yes','No','Yes');
spontaneousWindow = -flip(responseWindow);
if strcmpi(sponAns,'No')
    spontPrompt = "Time before the trigger in [s] (e.g. -0.8, -0.6 s)";
    sponDef = string(sprintf('%.3f, %.3f',spontaneousWindow(1),...
        spontaneousWindow(2)));
    sponStr = inputdlg(spontPrompt, 'Inputs',[1,30],sponDef);
    if ~isempty(sponStr)
        spontAux = str2num(sponStr{1});
        if length(spontAux) ~= 2 || spontAux(1) > spontAux(2) || ...
                spontAux(1) < timeLapse(1)
            fprintf(1, 'The given input was not valid.\n')
            fprintf(1, 'Keeping the mirror version!\n')
        else
            spontaneousWindow = spontAux;
        end
    end
end
fprintf(1,'Spontaneous window: %.2f to %.2f ms before the trigger\n',...
    spontaneousWindow(1)*1e3, spontaneousWindow(2)*1e3)


%% We now need to consider one stack of triggers that we will later split up into their respective conditions.


% Some conditions are duplicated as either '_AllTriggers', or 'Block'. The
% Block triggers are 5 seconds time durations of the 'All_Triggers'
% condition. 

    % e.g. 10Hz_All_Triggers will contain trigger times for each TTL pulse,
    % which can be a little as 5ms. The Block variety of this condition
    % considers the duration of the pulse train (usually 5 seconds), to compare
    % with other longer manipulations.

% We suggest going for a 'Mech_All_Block' conditions initally.

condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
condGuess = contains(condNames, 'whiskerall', 'IgnoreCase', true);


% Choose the conditions to create the stack upon

[chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
    'PromptString',...
    'Choose the condition which has all triggers: (one condition)',...
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

%% We'll now create 3-D stacks of the chosen data that we can access in alter analysis.


% discStack - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.

[discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);






% Ne = Number of units + the Mech as the first event + the Laser as the last
% event, 

% Nt = Number of time samples in between the time window, 

% Nta = number of total triggers.

[Ne, Nt, NTa] = size(discStack);




% Computing the time axis for the stack

tx = (0:Nt - 1)/fs + timeLapse(1);


%% Now we can choose the conditions we wish to analyse within our 'All' condition. 

% Choose the conditions to look at
auxSubs = setdiff(1:numel(condNames), chCond);
ccondNames = condNames(auxSubs);
[cchCond, iOk] = listdlg('ListString',ccondNames,'SelectionMode','multiple',...
    'PromptString',...
    'Choose the condition(s) to look at (should be in the a subset of ALL condition):',...
    'ListSize', [350, numel(condNames)*16]);
if ~iOk
    fprintf(1,'Cancelling...\n')
    return
end

 cchCond = flip(cchCond); % Flip the conditions to get Control on x-axis later

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

% Adding all the triggers from the piezo and the laser in one array
% allWhiskersPlusLaserControl = ...
%     union(Conditions(allWhiskerStimulus).Triggers,...
%     Conditions(laserControl).Triggers,'rows');


%% We'll create logical Boolean flags to index the subset of conditinos into the larger stack
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
[Results, Counts] = statTests(discStack, delayFlags, timeFlags);

indCondSubs = cumsum(Nccond:-1:1);
consCondNames = condNames(consideredConditions);
% Plotting statistical tests
[Figs, Results] = scatterSignificance(Results, Counts,...
    consCondNames, delta_t, sortedData(goods,1));
configureFigureToPDF(Figs);
stFigBasename = fullfile(figureDir,[expName,' ']);
stFigSubfix = sprintf(' Stat RW%.1f-%.1fms SW%.1f-%.1fms',...
    responseWindow(1)*1e3, responseWindow(2)*1e3, spontaneousWindow(1)*1e3,...
    spontaneousWindow(2)*1e3);
ccn = 1;
%for cc = indCondSubs
for cc = 1:numel(Figs)
    if ~ismember(cc, indCondSubs)
        altCondNames = strsplit(Figs(cc).Children(2).Title.String,': ');
        altCondNames = altCondNames{2};
    else
        altCondNames = consCondNames{ccn};
        ccn = ccn + 1;
    end
    stFigName = [stFigBasename, altCondNames, stFigSubfix];
    %     if ~exist([stFigName,'.pdf'],'file') || ~exist([stFigName,'.emf'],'file')
    %         print(Figs(cc),[stFigName,'.pdf'],'-dpdf','-fillpage')
    %         print(Figs(cc),[stFigName,'.emf'],'-dmeta')
    %     end
    savefig(Figs(cc),fullfile(figureDir, [altCondNames, '.fig']));
end

H = cell2mat(cellfun(@(x) x.Pvalues,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), 'UniformOutput', 0),...
    'UniformOutput', 0)) < 0.05;

Htc = sum(H,2);
CtrlCond = contains(consCondNames,'control','IgnoreCase',true);
if ~nnz(CtrlCond)
    CtrlCond = true(size(H,2),1);
end
wruIdx = any(H(:,CtrlCond),2);
Nwru = nnz(wruIdx);

fprintf('%d responding units:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})



%% Adding our newly calculated results to the clInfo table.


window = diff(responseWindow);

for a = 1: length(consCondNames)
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_Rate_Spont']} = mean(Counts{a,1}')'/window;
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_Rate_Evoked']} = mean(Counts{a,2}')'/window;
end


% Significant mechanical responses per condition
b = length(consCondNames);
for a = 1 : length(consCondNames)
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_R']} = Results(b).Activity(1).Pvalues < 0.05;
    b = b + length(consCondNames) - a;
end


sZ = size(Counts);
d = length(consCondNames);
e = 1;
f = length(consCondNames);
for a = 1:length(consCondNames) - 1
    for b = (a + 1): length(consCondNames)

        for c = 1: sZ(1,2)
            clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a},'_vs_',consCondNames{1,b}, '_', Results(e).Activity(c).Type, '_Response']} = Results(e).Activity(c).Pvalues < 0.05; 

        end

        e = e + 1;

        if e == f
            e = e + 1;
        end
    end
    d = d - 1;
    f = f + d;
end

%% Adding the absolute depth of each unit to the clInfo table.

dpth = 1500;

abs_depth = dpth - clInfo.depth;

clInfo = addvars(clInfo,abs_depth,'After','depth',...
    'NewVariableNames','abs_depth');


%% Opto-tagging

% Now that we have the absolute depths of our units, we can look at the
% spike latencies of each unit in response to the laser laser triggers, to
% determine the temporal sequence of spike activation (e.g. which
% population is expressing the channel-rhodopsin protein).

[TaggedIDs, fig] = Optotag([2,10], [0.2,4], clInfo, gclID, Conditions(14).name, Conditions(14).Triggers, sortedData, fs);

idxTagged = ismember(clInfo.id, TaggedIDs);
if ~any(ismember(clInfo.Properties.VariableNames,'Tagged'))
    clInfo = addvars(clInfo,idxTagged,'After','ActiveUnit',...
        'NewVariableNames','Tagged');
end

%% Saving the ammended clInfo table

writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));


%% Adding the mean trace signals to the Conditions variable
if ~isfield(Conditions,'Stimulus') ||...
        any(arrayfun(@(x) isempty(x.Stimulus), Conditions(consideredConditions)))
    fprintf(1,'Writting the stimulus raw signal into Conditions variable:\n')
    whFlag = contains(trigNames, whStim, 'IgnoreCase', 1);
    lrFlag = contains(trigNames, cxStim, 'IgnoreCase', 1);
    cdel = 1;
    for cc = consideredConditions
        fprintf(1,'- %s\n', Conditions(cc).name)
        Conditions(cc).Stimulus = struct(...
            'Mechanical',reshape(mean(cst((1),:,delayFlags(:,cdel)),3),...
            1,Nt),...
            'MechPressure',reshape(mean(cst((2),:,delayFlags(:,cdel)),3),...
            1,Nt),...
            'Laser',reshape(mean(cst(lrFlag,:,delayFlags(:,cdel)),3),...
            1,Nt),'TimeAxis',(0:Nt-1)/fs + timeLapse(1));
        cdel = cdel + 1;
    end
end


%% Creating a configuration structure to keep track of all the 
configStructure = struct('Experiment', fullfile(dataDir,expName),...
    'Viewing_window_s', timeLapse, 'Response_window_s', responseWindow,...
    'BinSize_s', binSz, 'Trigger', struct('Name', condNames{chCond},...
    'Edge',onOffStr), 'ConsideredConditions',{consCondNames});


%% Filter question
filterIdx = true(Ne,1);
ansFilt = questdlg('Would you like to filter for significance?','Filter',...
    'Yes','No','Yes');
filtStr = 'unfiltered';
if strcmp(ansFilt,'Yes')
    filterIdx = [true; wruIdx];
    filtStr = 'Filtered';
end
%ruIdx = wruIdx;


%% Getting the relative spike times for the responsive units (ru)
% For each condition, the first spike of each ru will be used to compute
% the standard deviation of it.
cellLogicalIndexing = @(x,idx) x(idx);
isWithinResponsiveWindow =...
    @(x) x > responseWindow(1) & x < responseWindow(2);

firstSpike = zeros(Nwru,Nccond);
M = 16;
binAx = responseWindow(1):binSz:responseWindow(2);
condHist = zeros(size(binAx,2)-1, Nccond);
firstOrdStats = zeros(2,Nccond);
condParams = zeros(M,3,Nccond);
txpdf = responseWindow(1):1/fs:responseWindow(2);
condPDF = zeros(numel(txpdf),Nccond);
csvBase = fullfile(dataDir, expName);
csvSubfx = sprintf(' VW%.1f-%.1f ms.csv', timeLapse(1)*1e3, timeLapse(2)*1e3);
existFlag = false;
condRelativeSpkTms = cell(Nccond,1);
relativeSpkTmsStruct = struct('name',{},'SpikeTimes',{});
spkDir = fullfile(dataDir, 'SpikeTimes');
for ccond = 1:size(delayFlags,2)
    csvFileName = [csvBase,' ',consCondNames{ccond}, csvSubfx];
    relativeSpikeTimes = getRasterFromStack(discStack,~delayFlags(:,ccond),...
        filterIdx(3:end), timeLapse, fs, true, false);
    relativeSpikeTimes(:,~delayFlags(:,ccond)) = [];
    relativeSpikeTimes(~filterIdx(2),:) = [];
    condRelativeSpkTms{ccond} = relativeSpikeTimes;

    clSpkTms = cell(size(relativeSpikeTimes,1),1);
    if exist(csvFileName, 'file') && ccond == 1
        existFlag = true;
        ansOW = questdlg(['The exported .csv files exist! ',...
            'Would you like to overwrite them?'],'Overwrite?','Yes','No','No');
        if strcmp(ansOW,'Yes')
            existFlag = false;
            fprintf(1,'Overwriting... ');
        end
    end
    fID = 1;

    for cr = 1:size(relativeSpikeTimes, 1)
        clSpkTms(cr) = {sort(cell2mat(relativeSpikeTimes(cr,:)))};
        if fID > 2
            fprintf(fID,'%s,',gclID{cr});
            fprintf(fID,'%f,',clSpkTms{cr});fprintf(fID,'\n');
        end
    end
    if fID > 2
        fclose(fID);
    end
    relativeSpkTmsStruct(ccond).name = consCondNames{ccond};
    relativeSpkTmsStruct(ccond).SpikeTimes = condRelativeSpkTms{ccond};
  
end
save(fullfile(dataDir,[expName,'_exportSpkTms.mat']),...
    'relativeSpkTmsStruct','configStructure')


%% Ordering PSTH

orderedStr = 'ID ordered';
dans = questdlg('Do you want to order the PSTH other than by IDs?',...
    'Order', 'Yes', 'No', 'No');
ordSubs = 1:nnz(filterIdx(2:Ncl+1));
pclID = gclID(filterIdx(2:end));
if strcmp(dans, 'Yes')
    if ~exist('clInfo','var')
        clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
    end
    % varClass = varfun(@class,clInfo,'OutputFormat','cell');
    [ordSel, iOk] = listdlg('ListString', clInfo.Properties.VariableNames,...
        'SelectionMode', 'multiple');
    orderedStr = [];
    ordVar = clInfo.Properties.VariableNames(ordSel);
    for cvar = 1:numel(ordVar)
        orderedStr = [orderedStr, sprintf('%s ',ordVar{cvar})]; %#ok<AGROW>
    end
    orderedStr = [orderedStr, 'ordered'];
    
    if ~strcmp(ordVar,'id')
        [~,ordSubs] = sortrows(clInfo(pclID,:),ordVar, 'ascend');
    end
end

%% Plot PSTH
goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
psthFigs = gobjects(Nccond,1);
for ccond = 1:Nccond
    figFileName = sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, Conditions(consideredConditions(ccond)).name, timeLapse*1e3,...
        binSz*1e3, responseWindow*1e3, spontaneousWindow*1e3, onOffStr,...
        orderedStr, filtStr);
    [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx,:,:),timeLapse,...
        ~delayFlags(:,ccond),binSz,fs);
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    psthFigs(ccond) = plotClusterReactivity(PSTH(ordSubs,:,ccond),trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};...
        pclID(ordSubs)],...
        strrep(expName,'_','\_'),...
        stims, csNames);
    configureFigureToPDF(psthFigs(ccond));
    psthFigs(ccond).Children(end).YLabel.String =...
        [psthFigs(ccond).Children(end).YLabel.String,...
        sprintf('^{%s}',orderedStr)];

    
end
% for a = 1:length(consideredConditions)
%     savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH_0.001binSz.fig']));
% end


%% Plot depthPSTH

csNames = fieldnames(Triggers);
pclTblInd = ismember(clInfo.id, pclID(ordSubs));
depthOrd = sort(clInfo.abs_depth(pclTblInd), 'descend');
depthVals = [4200:-200:0];
depthInds = zeros(length(depthVals), 1);
for a = 1:length(depthInds)
    selDepths = depthOrd < depthVals(a);
    if sum(selDepths)> 0
        depthInd = find(selDepths);
        depthInds(a) = depthInd(1);
    end
end
topDepth = find(depthInds == 0);
topDepth = topDepth(1);
depthInds(topDepth) = length(pclID);
depthInds = depthInds(depthInds > 0);
selClIDs = cell(length(pclID),1);
depthVals = depthVals(1:length(depthInds))';
for a = 1:length(depthInds)
    selClIDs{depthInds(a)} = num2str(depthVals(a));
end
goodsIdx = logical(clInfo.ActiveUnit);
%csNames = csNames(2:end);
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
psthFigs = gobjects(Nccond,1);
for ccond = 1:Nccond
    figFileName = sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, Conditions(consideredConditions(ccond)).name, timeLapse*1e3,...
        binSz*1e3, responseWindow*1e3, spontaneousWindow*1e3, onOffStr,...
        orderedStr, filtStr);
    [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx,:,:),timeLapse,...
        ~delayFlags(:,ccond),binSz,fs);
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    %stims = stims([2,3],:);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    psthFigs(ccond) = plotClusterReactivityRoss(PSTH(ordSubs,:,ccond),trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};...
        selClIDs],...
        strrep(expName,'_','\_'),...
        stims, csNames);
    configureFigureToPDF(psthFigs(ccond));
    psthFigs(ccond).Children(end).Title.String = consCondNames{ccond};
   
end
% for a = 1:length(consideredConditions)
%     savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH_0.001binSz.fig']));
% end



%% We can now use k means clustering analysis for a less biased splitting of the population based on the units' depths, latencies and jitter.
[GroupIDs, fig] = findSubPopulations(3, clInfo, gclID, Conditions(14).name, Conditions(14).Triggers, sortedData, fs);


Group1 = GroupIDs{1,2};
Group2 = GroupIDs{2,2};
Group3 = GroupIDs{3,2};

idxGroup1 = ismember(gclID, Group1);
idxGroup2 = ismember(gclID, Group2);
idxGroup3 = ismember(gclID, Group3);

% Order the groups by depth

ruIdx = [true, true, true; idxGroup3, idxGroup2, idxGroup1];
%% We can plot the Depth vs change in FR by the groups we've just assigned


% Prep Variables
mrControl = clInfo.Mech_Control_4mW_Mech_Rate_Evoked-clInfo.Mech_Control_4mW_Mech_Rate_Spont;
mrLaser = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked-clInfo.Mech_Laser_5sec_4mW_Mech_Rate_Spont;
data = table(clInfo.id, clInfo.abs_depth, mrControl, mrLaser);
mrInd = ismember(data.Var1, clInfo.id(clInfo.Mech_Control_4mW_R==true));
data = data(find(mrInd),:);
data = sortrows(data,'Var2','ascend');
ID = data.Var1;
depths = data.Var2;
c = data.Var3;
l = data.Var4;
rng('default');
jitter = randi(20,size(depths));
depths = -1*(depths + jitter);

red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];
purple = [0.5,0,0.5];
grey = [0.5, 0.5, 0.5];
colour = [green; red; blue];

% Depth vs Evoked Rate (Control): MR Evoked - MR Baseline


deltaRate = zeros(length(depths), 2);
deltaRate(:,2) = c;

figure('Color', 'White', 'Name', 'Changes in FR with Mechanical Stimulation');
hold on
for unit = 1:length(ID)
dr = deltaRate(unit,:);
dpth = [depths(unit), depths(unit)];
if ismember(ID{unit}, Group3)
plot(dr, dpth, 'LineStyle', '-', 'Color', [green, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Group2)
plot(dr, dpth, 'LineStyle', '-', 'Color', [red, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Group1)
plot(dr, dpth, 'LineStyle', '-', 'Color', [blue, 0.5], 'LineWidth', 10);
end
end
ax = gca;
ax.Title.String = 'Changes in FR with Mechanical Stimulation';
ax.YLabel.String = 'Depth [\mum]';
ax.XLabel.String = '\DeltaRate [Hz]';
ax.FontName = 'Arial';
ax.FontSize = 30;
ax.YLim = [-1600, 0];
% ax.XLim = [-50, 50];
% ax.XTick = [-50:10:50];
plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');


% Mech vs Mech Laser: MR Laser - MR Evoked


deltaRate = zeros(length(depths), 2);
deltaRate(:,2) = l-c;

figure('Color', 'White', 'Name', 'Laser Modulation of Mechanical Responses');
hold on
for unit = 1:length(ID)
dr = deltaRate(unit,:);
dpth = [depths(unit), depths(unit)];
if ismember(ID{unit}, Group3)
plot(dr, dpth, 'LineStyle', '-', 'Color', [green, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Group2)
plot(dr, dpth, 'LineStyle', '-', 'Color', [red, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Group1)
plot(dr, dpth, 'LineStyle', '-', 'Color', [blue, 0.5], 'LineWidth', 10);
end
end
ax = gca;
ax.Title.String = 'Laser Modulation of Mechanical Responses';
ax.YLabel.String = 'Depth [\mum]';
ax.XLabel.String = '\DeltaRate [Hz]';
ax.FontName = 'Arial';
ax.FontSize = 30;
ax.YLim = [-1600, 0];
% ax.XLim = [-50, 50];
% ax.XTick = [-50:10:50];
plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');




%% We can also plot the popPSTHs according to these groups 

red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];
purple = [0.5,0,0.5];
grey = [0.5, 0.5, 0.5];
csNames = fieldnames(Triggers);
colour = [green; red; blue];
% ruIdx = filterIdx;
axLabel = 'Firing Rate [Hz]';
for ccond = 1: length(consideredConditions)
    ruIdxDim = size(ruIdx);
    nFilters = ruIdxDim(2);
    
    fig = figure('Name',[expName,'_', consCondNames{ccond}],'Color',[1,1,1]);
    fthAxFlag = false;
    if ~exist('stims','var')
        totlX = nFilters;
        
    elseif ~isempty(stims)
        fthAxFlag = true;
        totlX = nFilters + 1;
    end
   
    
    clrDivider = nFilters-1;
    if nFilters <= 1
        clrDivider = nFilters;
    end
   
    datclr = [green; red; blue];   
    
    for a = 1:nFilters
        ax(a) = subplot(totlX,1,a,'Parent',fig);
        [PSTHarray{a}(:,:,ccond), trig, sweeps] = getPSTH(discStack(ruIdx(:,a),:,:),timeLapse,...
            ~delayFlags(:,ccond),binSz,fs);
        
        PSTH = PSTHarray{a}(:,:,ccond);
        [Ncl, Npt] = size(PSTH);
        PSTHn = PSTH./max(PSTH,[],2);
        
        
        psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
        trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
        popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps * binSz);
        plot(psthTX,popPSTH,'Color', colour(a,:))
        ax(a).XLabel.String = sprintf('Time_{%.2f ms} [s]',binSz*1e3);
        ax(a).XLim = [timeLapse(1), timeLapse(2)];

        ax(a).Box = 'off';
        ax(a).ClippingStyle = 'rectangle';
        legend(ax(a),'show','Location','northeast')
        ruNames = [{'Group 3'}, {'Group 2'}, {'Group 1'}];
        ax(a).Legend.String = ruNames{a};
        ax(a).YAxis(1).Label.String = axLabel;
        
    end
    ax(1).Title.String = consCondNames(ccond);

    ax2 = subplot(totlX,1,nFilters + 1,'Parent',fig);
    
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    
    [r,c] = size(stims);
    
    if r < c
        stims = stims';
        stmClr = zeros(r, 3);
        
    end
    
    
    stmClr([1,5,9]) = 1;
    stmClr(stmClr == 0) = 0.25;
    
    
    
    for cs = 1:min(r,c)
        if exist('IDs','var')
            plot(ax2,trigTX,stims(:,cs),'LineStyle','-.','LineWidth',0.5,...
                'DisplayName', csNames{cs}, 'Color', stmClr(cs,:))
        else
            plot(ax2,trigTX,stims(:,cs),'LineStyle','-.','LineWidth',0.5,  'Color', stmClr(cs,:))
        end
        
        
        if cs == 1
            ax2.NextPlot = 'add';
        end
    end
    legend(ax2,'show','Location','best')
    
    ax2.Box = 'off';
    ax2.XLim = [timeLapse(1), timeLapse(2)];
    ax2.XAxis.Visible = 'off';
    ax2.YAxis.Visible = 'off';
    linkaxes([ax, ax2], 'x')
    %end
    
    
 
    %
    % Formatting the population PSTH plot
    
    ax2.YAxis(1).Limits = [0,1.01];
    ax2.Legend.String = csNames;
    %
    %
    
    %
    %
    clear PSTH
    clear PSTHarray
    clear ax
    
end
% end

clear ax
%% We can now look at individual units and their spiking responses to each stimulus over n trials in a Spike Raster plot

csNames = fieldnames(Triggers);
% csNames = csNames(2:end);
IDs = csNames;
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));

rasterDir = fullfile(figureDir,'Rasters\');
if ~mkdir(rasterDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end


Power = NaN(length(consCondNames),1);
for cc = 1:length(consCondNames)
    
    mWfind = strfind(consCondNames{cc}, 'mW');
    
    pwr = consCondNames{cc}(mWfind-2:mWfind+1);
    if isnan(pwr)
        pwrMissing = true;
    elseif contains(pwr, '.')
        pwr = consCondNames{cc}(mWfind-3:mWfind+1);
    elseif contains(pwr, '_') || contains(pwr, ' ')
        pwr = pwr(2:end);
    end
    pwr = str2double(pwr(1:end-2));
    Power(cc) = pwr;
end
pwrs = unique(Power);



for a = 1:length(pwrs)
    pwr = pwrs(a);

    pwrInd = Power == pwr;
    clIDind =  clInfo.id(clInfo.Mech_Control_4mW_R==true);
    lngth = length(clIDind);
    for a = 1:lngth
        rng('default');
        cl = clIDind(a);
        clSel = find(ismember(pclID, cl));

        rasCondSel = flip(find(pwrInd));
        rasCond = consideredConditions(rasCondSel);
        rasCondNames = consCondNames(rasCondSel);
        Nrcl = numel(clSel);
        % Reorganize the rasters in the required order.
        clSub = find(ismember(gclID, pclID(clSel)))+1;
        [rasIdx, rasOrd] = ismember(pclID(ordSubs), pclID(clSel));
        clSub = clSub(rasOrd(rasIdx));
        clSel = clSel(rasOrd(rasOrd ~= 0));
        Nma = min(Na(rasCondSel));
        rasFig = figure;
        columns = length(pwrs);
        Nrcond = length(rasCond);
        ax = gobjects(Nrcond*Nrcl,1);
        for cc = 1:length(rasCond)
            % Equalize trial number
            trigSubset = sort(randsample(Na(rasCondSel(cc)),Nma));
            tLoc = find(delayFlags(:,rasCondSel(cc)));
            tSubs = tLoc(trigSubset);
            % Trigger subset for stimulation shading
            trigAlSubs = Conditions(rasCond(cc)).Triggers(trigSubset,:);
            timeDur = round(diff(trigAlSubs, 1, 2)/fs, 3);
            trigChange = find(diff(timeDur) ~= 0);
            
            for ccl = 1:Nrcl
                lidx = ccl + (cc - 1) * Nrcl;
                ax(lidx) = subplot(Nrcond, Nrcl, lidx);       %  subplot( Nrcl, Nrcond, lidx); % to plot the other way around
                title(ax(lidx),sprintf(rasCondNames{cc}), 'Interpreter', 'none') % ,pclID{clSel(ccl)}
                plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                    timeLapse, fs,'',ax(lidx));
                ax(lidx).YAxisLocation = 'origin';ax(lidx).YAxis.TickValues = Nma;
                ax(lidx).YAxis.Label.String = 'Trials';

                xlabel(ax(lidx), 'Time [s]')
                initSub = 0;
                optsRect = {'EdgeColor','none','FaceColor','none'};
                for ctr = 1:numel(trigChange)
                    rectangle('Position',[0, initSub,...
                        timeDur(trigChange(ctr)), trigChange(ctr)],optsRect{:})
                    initSub = trigChange(ctr);
                end
                rectangle('Position', [0, initSub, timeDur(Nma),...
                    Nma - initSub],optsRect{:})
                
                
                stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
%                 stims = stims([2,3],:);
                
                stims = stims - median(stims,2);
                
                
                
                for cs = 1:size(stims,1)
                    if abs(log10(var(stims(cs,:),[],2))) < 13
                        [m,b] = lineariz(stims(cs,:),1,0);
                        stims(cs,:) = m*stims(cs,:) + b;
                    else
                        stims(cs,:) = zeros(1,Nt);
                    end
                end
                
                [r,c] = size(stims);
                
                if r < c
                    stims = stims';
                    stmClr = zeros(r, 3);
                    
                end
                
                hold off
                ax = gca;
                
               

                ax.FontName ='Arial';
                ax.FontSize = 25;
                      
                
                
            end
        end
        
        
        
        
        
        rasConds = rasCondNames{1};
        if length(rasCondNames) > 1
            for r = 2:length(rasCondNames)
                rasConds = [rasConds, '+', rasCondNames{r}];
            end
        end
        
        linkaxes(ax,'x')
        rasFigName = ['Unit_', cell2mat(cl), '_', ];
        rasFig.Name = [rasFigName, '_', num2str(pwr), 'mW'];
        configureFigureToPDF (rasFig);
        set(rasFig, 'Position', get(0, 'ScreenSize'));
        saveas(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds,'_', num2str(timeLapse(1)), '_to_', num2str(timeLapse(2)),'.emf']));
        %savefig(rasFig,fullfile(rasterDir, [rasFigName, ' ', num2str(pwr), 'mW.fig']));
        savefig(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds, '.fig']));
    end
end
% close all

%% ISI analysis - 


% Computing the Triggered ISIs


% Getting the ActiveUnit ISIs from sortedData
ind = ismember(gclID,(clInfo.id(clInfo.Mech_Control_4mW_R==true)));
% gclID = clInfo.id(ind == true);
% ind = ismember(sortedData(:,1), gclID(wruIdx,1));
% ind = ismember(sortedData(:,1), gclID(wru));
Ncl = sum(ind);
spkSubs = cellfun(@(x) round(x.*fs), sortedData(ind,2),...
    'UniformOutput', false);
ISIVals = cellfun(@(x) [x(1)/fs; diff(x)/fs], spkSubs,...
    'UniformOutput', 0);
NnzvPcl = cellfun(@numel,ISIVals);
Nnzv = sum(NnzvPcl);
rows = cell2mat(arrayfun(@(x,y) repmat(x,y,1), (1:Ncl)', NnzvPcl,...
    'UniformOutput', 0));
cols = cell2mat(spkSubs);
vals = cell2mat(ISIVals);
ISIspar = sparse(rows, cols, vals);

% Creating the TrigISIs Struct


ConsConds = cchCond;
nCond = length(ConsConds);
% sortedData = sortedData(:,1);
for chCond = 1:nCond
    TrigISIs.name = Conditions(cchCond(chCond)).name;
    TrigISIs.Vals(1).name = 'Spontaneous';
    TrigISIs.Vals(2).name = 'Evoked';


% Adding ISIs to TrigISIs
spontaneousWindow = -flip(responseWindow);


    
    
    % contains will give multiple units when looking for e.g. cl45
%     spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % spkSubs replaces round(sortedData{goods(1),2}*fs) for the rest of the
    % clusters
    % Subscript column vectors for the rest good clusters
    for wIndex = 1:2
        if wIndex == 1
            Window = spontaneousWindow;
        else
            Window = responseWindow;
        end
        [~, isiStack] = getStacks(false,Conditions(cchCond(chCond)).Triggers, onOffStr,...
            Window,fs,fs,[],ISIspar);
        lInda = isiStack > 0; 
        % timelapse becomes spontaneousWindow for pre-trigger, and responseWindow
        % for post
        TrigISIs.Vals(wIndex).TriggeredIsI = isiStack;
        for histInd = 1: Ncl
            figure('Visible','off');
            hisi = histogram(log10(isiStack(histInd,:,:)), 'BinEdges', log10(0.0001):0.01:log10(10));
            TrigISIs.Vals(wIndex).cts{histInd} = hisi.BinCounts;
            TrigISIs.Vals(wIndex).bns{histInd} = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
            
            close gcf;
        end
    end

% ISIs and CumISIs

    for a = 1:length(TrigISIs.Vals(wIndex).cts)
        TrigISIs.Vals(1).ISI{a} = TrigISIs.Vals(1).cts{a}./sum(TrigISIs.Vals(1).cts{a});
        TrigISIs.Vals(2).ISI{a} = TrigISIs.Vals(2).cts{a}./sum(TrigISIs.Vals(2).cts{a});
        TrigISIs.Vals(1).CumISI{a} = cumsum(TrigISIs.Vals(1).ISI{a});
        TrigISIs.Vals(2).CumISI{a} = cumsum(TrigISIs.Vals(2).ISI{a});
    end



% Plotting the ISIs


figure('Color', 'White', 'Name', [TrigISIs.name]);
hold on
for wIndex = 1:2
    CondNames = {'Baseline', 'Evoked'};
    clr = {[0.5, 0, 0], [0, 1, 1]};
    IsiStack = zeros(length(TrigISIs.Vals(wIndex).ISI), length(TrigISIs.Vals(wIndex).ISI{1}));
    cumIsiStack = zeros(length(TrigISIs.Vals(wIndex).CumISI), length(TrigISIs.Vals(wIndex).CumISI{1}));
    for cInd = 1:length(TrigISIs.Vals(wIndex).CumISI)
        IsiStack(cInd,:) =  TrigISIs.Vals(wIndex).ISI{cInd};
        cumIsiStack(cInd,:) =  TrigISIs.Vals(wIndex).CumISI{cInd};
    end
    % Getting rid of NaNs
    cumIsiStack(length(TrigISIs.Vals(wIndex).CumISI) + 1,:) = NaN;
    IsiStack(length(TrigISIs.Vals(wIndex).ISI) + 1,:) = NaN;
    cumIndNaN = (length(TrigISIs.Vals(wIndex).CumISI) + 1);
    for nInd = 1:(length(TrigISIs.Vals(wIndex).CumISI))
        if isnan(cumIsiStack(nInd,:)) == true
            cumIndNaN = [cumIndNaN; nInd];
        end
    end
    cumIndNaN = sort(cumIndNaN, 'descend');
    cumIsiStack(cumIndNaN,:) =[];
    IsiStack(cumIndNaN,:) = [];
    stckSz = size(cumIsiStack);
    ISIcum = sum(cumIsiStack)/stckSz(1);
    % Only if the bin widths are constant!!!
    plot(TrigISIs(1).Vals(1).bns{1}, ISIcum, 'Color', clr{wIndex})
end
legend(CondNames)
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([-4, 1]);
xticks([-4:1])
ylim([0, 1]);
fig = gcf;
ax = gca;
ax.FontSize = 20;
ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;
end
%%
% We've now come to the end of the tutorial script. Please feel free to try
% uot other conditions, bin sizes, optotagging paramters etc. Just play
% around with it until you get bored.
