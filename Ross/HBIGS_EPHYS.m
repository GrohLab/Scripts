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

% probe file = 'Z:\Ross\ProbeFiles\Corrected_H3_ChanMap.mat'
% sampling frequency = 3.003003003003003e+04




%% Once sorted and curated we can import our units and their spikes into MATLAB.

cd C:\Users\NeuroNetz\Documents\GitHub\Scripts\Ross

% Experiment Name
expName = 'Replace_this_text_with_your_experiment_name';

% Spike-sorted units and their spike times
sortedData = importPhyFiles(dataDir, expName);

%Unit/Cluster (same thing) information in a table.
clInfo = getClusterInfo([dataDir, '\cluster_info.tsv']);


%% Now we need to organise the TTL data to know when we are stimulating the brain.

% Calling a function to give us the Conditions and the Trigger times.
[Conditions, Triggers] = getTriggers(expName, dataDir, fs);


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
defInputs = {'-2, 6', '0.1, 1.5', '0.01'};
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


% Some conditinon are duplicated as eithe '_AllTriggers', or 'Block'. The
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

[Tagged, fig] = Optotag([2,10], [0.2,4], clInfo, gclID, Conditions(21).name, Conditions(21).Triggers, sortedData, fs);

idxTagged = ismember(clInfo.id, TaggedIDs);
if ~any(ismember(clInfo.Properties.VariableNames,'Tagged'))
    clInfo = addvars(clInfo,idxTagged,'After','ActiveUnit',...
        'NewVariableNames','Tagged');
end

%% Saving the ammended clInfo table

writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));


%% Addition mean signals to the Conditions variable
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


%% Configuration structure
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