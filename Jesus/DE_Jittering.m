% 30.08.19 Jittering analysis by using the Data Explorer.
clearvars
%% Selecting experiment
% Choosing the working directory
dataDir = uigetdir('E:\Data\VPM\Jittering\Silicon Probes\',...
    'Choose a working directory');
if dataDir == 0
    return
end
%% Loading data
% Creating the figure directory
FigureDir = fullfile(dataDir,'Figures\');
if ~exist(FigureDir, "dir")
    if ~mkdir(FigureDir)
        error("Could not create figure directory!\n")
    end
end
%WARNING! Protocol getter is not fully stable!
if isempty(dir(fullfile(dataDir, '*analysis.mat')))
    pgObj = ProtocolGetter(dataDir);
    pgObj.getConditionSignals;
    pgObj.getSignalEdges;
    pgObj.getFrequencyEdges;
    pgObj.pairStimulus;
    pgObj.saveConditions;
end
% Loading the necessary files
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables\n')
    return
end
fnOpts = {'UniformOutput', false};
axOpts = {'Box','off','Color','none'};
spk_file_vars = {'spike_times','gclID','Nt','Ns','goods'};
%% Constructing the helper 'global' variables

spkPttrn = "%s_Spike_Times.mat";
spk_path = fullfile(dataDir, sprintf(spkPttrn, expName));
if ~exist(spk_path, "file")
    % Number of total samples
    Ns = structfun(@numel,Triggers); Ns = min(Ns(Ns>1));
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
    badsIdx = badsIdx | silentUnits;
    if ~any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
        clInfo = addvars(clInfo,~badsIdx,'After',1,...
            'NewVariableNames','ActiveUnit');
        writeClusterInfo(clInfo, fullfile(dataDir, 'cluster_info.tsv'), 1);
    end
    gclID = sortedData(goods,1);
    badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
    spike_times = sortedData(~badsIdx, 2);
    save(spk_path, spk_file_vars{:})
else
    % Subscript column vectors for the rest good clusters
    load(spk_path, spk_file_vars{:})
end
spkSubs = cellfun(@(x) round(x.*fs),spike_times, fnOpts{:});
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

%% Inter-spike intervals
isiFile = fullfile(dataDir,[expName,'_ISIvars.mat']);
if ~exist(isiFile,'file')
    spkSubs2 = cellfun(@(x) round(x.*fs), sortedData(goods,2),...
        fnOpts{:});
    ISIVals = cellfun(@(x) [x(1)/fs; diff(x)/fs], spkSubs2,...
        fnOpts{:});
    NnzvPcl = cellfun(@numel,ISIVals);
    Nnzv = sum(NnzvPcl);
    rows = cell2mat(arrayfun(@(x,y) repmat(x,y,1), (1:Ncl)', NnzvPcl,...
        fnOpts{:}));
    cols = cell2mat(spkSubs2);
    vals = cell2mat(ISIVals);
    try
        ISIspar = sparse(rows, cols, vals);
    catch
        fprintf(1, 'Not possible to create such a big array')
    end
else
    % load(isiFile,'ISIspar')
end
% ISIsignal = zeros(Ncl,Ns,'single');
% for ccl = 1:Ncl
%     ISIsignal(ccl,spkSubs2{ccl}) = ISIVals{ccl};
% end

%% User controlling variables
% Parameter configuration (pc) file name (FN), file path (FP) and variables
% to save (V2S)
pcPttrn = "%s_ParameterConfiguration.mat";
pcFN = sprintf(pcPttrn, expName); pcFP = fullfile(dataDir, pcFN);
pcV2S = {'configStructure', 'confKeys'}; pcSub = []; pcFlag = false;
if exist(pcFP, 'file')
    pcFlag = true;
    load(pcFP, pcV2S{:});
    pcSub = listdlg("ListString", join(confKeys,2), ...
        "PromptString", ...
        "Viewing       Response       Spontaneous     Trigger  "+ ...
        "Conditions      On/Off  BinSize [ms]", ...
        "name", "Configuration selection", ...
        "ListSize", [700, size(confKeys,1)*25], ...
        "CancelString", "New");
end

if isempty(pcSub)
    % If user chose to create a new parameter configuration set or the file
    % doesn't exist
    if ~pcFlag
        % File doesn't exist. Create it.
        configStructure = []; confKeys = [];
    else
        pcV2S = [pcV2S(:)', {'-append'}];
    end
    %% Select Time lapse, bin size, and spontaneous and response windows
    promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
        'Bin size [s]:'};
    defInputs = {'-0.1, 0.1', '0.002, 0.05', '0.001'};
    answ = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);
    if isempty(answ)
        fprintf(1,'Cancelling...\n')
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
    fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse*1e3)
    fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow*1e3)
    fprintf(1,'Bin size: %.3f ms\n', binSz*1e3)
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
    %% Select a condition to build the stack
    condNames = arrayfun(@(x) string(x.name),Conditions);
    condGuess = contains(condNames, ["whiskerall";"puffall"], ...
        'IgnoreCase', true);
    % Choose the conditions to create the stack upon
    [chCond, iOk] = listdlg('ListString',condNames,'SelectionMode','single',...
        'PromptString',...
        'Choose the condition which has all whisker triggers: (one condition)',...
        'InitialValue', find(condGuess), 'ListSize', [350, numel(condNames)*16]);
    if ~iOk
        fprintf(1,'Cancelling...\n')
        return
    end
    %% Select the onset or the offset of a trigger
    fprintf(1,'Condition ''%s''\n', Conditions(chCond).name)
    onOffStr = questdlg('Trigger on the onset or on the offset?','Onset/Offset',...
        'on','off','Cancel','on');
    if strcmpi(onOffStr,'Cancel')
        fprintf(1,'Cancelling...\n')
        return
    end
    %% Considered conditions selection
    % Choose the conditions to look at
    auxSubs = setdiff(1:numel(condNames), chCond);
    ccondNames = condNames(auxSubs);
    [cchCond, iOk] = listdlg('ListString',ccondNames, ...
        'SelectionMode', 'multiple', 'PromptString',...
        'Choose the condition(s) to look at (including whiskers):',...
        'ListSize', [350, numel(condNames)*16]);
    if ~iOk
        fprintf(1,'Cancelling...\n')
        return
    end
    consCondSubs = auxSubs(cchCond);
    consCondNames = condNames(consCondSubs);
    %% Keys creation for parameter configuration comparison
    CC_key = '';
    for cc = 1:length(consCondNames)-1
        CC_key = [CC_key, sprintf('%s-', consCondNames{cc})]; %#ok<AGROW>
    end
    CC_key = [CC_key, sprintf('%s', consCondNames{end})];
    C_key = condNames{chCond};
    SW_key = sprintf('SW%.2f-%.2f', spontaneousWindow*1e3);
    RW_key = sprintf('RW%.2f-%.2f', responseWindow*1e3);
    O_key = sprintf('%s', onOffStr);
    currentMapKey = {RW_key, SW_key, C_key, CC_key, O_key};
    VW_key = sprintf("VW%.2f-%.2f", timeLapse*1e3);
    BZ_key = sprintf("BZ%.2f",binSz*1e3);
    %% Appending the new pc to the array and saving
    configStructure = [configStructure, struct('Experiment', ...
        fullfile(dataDir,expName), 'Viewing_window_s', timeLapse, ...
        'Response_window_s', responseWindow, ...
        'Spontaneous_window_s', spontaneousWindow, 'BinSize_s', binSz, ...
        'Trigger', struct('Name', condNames{chCond}, 'Edge',onOffStr), ...
        'ConsideredConditions',{consCondNames})];
    newConfKey = [VW_key, string(currentMapKey), BZ_key];
    confKeys = [confKeys; newConfKey];
    save(fullfile(dataDir, pcFN), pcV2S{:})
    configStructure = configStructure(end);
else
    %% Loading paramters to the workspace
    % User chose a previous configuration
    configStructure = configStructure(pcSub);
    condNames = arrayfun(@(x) string(x.name), Conditions);
    timeLapse = configStructure.Viewing_window_s;
    responseWindow = configStructure.Response_window_s;
    spontaneousWindow = configStructure.Spontaneous_window_s;
    binSz = configStructure.BinSize_s;
    chCond = find(condNames == string(configStructure.Trigger.Name));
    onOffStr = configStructure.Trigger.Edge;
    [~, consCondSubs] = find(condNames == ...
        string(configStructure.ConsideredConditions)');
    consCondNames = condNames(consCondSubs);
    clearvars pcSub
    %% Keys creation for parameter configuration comparison
    CC_key = '';
    for cc = 1:length(consCondNames)-1
        CC_key = [CC_key, sprintf('%s-', consCondNames{cc})]; %#ok<AGROW>
    end
    CC_key = [CC_key, sprintf('%s', consCondNames{end})];
    C_key = condNames{chCond};
    SW_key = sprintf('SW%.2f-%.2f', spontaneousWindow*1e3);
    RW_key = sprintf('RW%.2f-%.2f', responseWindow*1e3);
    O_key = sprintf('%s', onOffStr);
    currentMapKey = {RW_key, SW_key, C_key, CC_key, O_key};
    VW_key = sprintf("VW%.2f-%.2f", timeLapse*1e3);
    BZ_key = sprintf("BZ%.2f",binSz*1e3);
end
%% Creating ephys figure folder
subFigDir = sprintf("Ephys %s %s %s", VW_key, RW_key, SW_key);
subFigDir = fullfile(FigureDir, subFigDir);
ephFigDir = subFigDir;
metaNameFlag = false;
if ~exist(subFigDir, "dir")
    if ~mkdir(subFigDir)
        fprintf(1, "Couldn't create %s!\n", subFigDir)
        fprintf(1, "Keeping metadata in figure file names.\n")
        metaNameFlag = true;
        ephFigDir = FigureDir;
    end
end
%% Constructing the stack out of the user's choice
% discStack - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.
[discStack, cst] = getStacks(false,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);
discStack(2,:,:) = [];
% ISI stack
try
    [~, isiStack] = getStacks(spkLog,Conditions(chCond).Triggers, onOffStr,...
        timeLapse,fs,fs,[],ISIspar);
catch
    fprintf(1,'Not able to do the ISI stack\n')
end
% [dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
%     'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
%     continuousSignals);
if ~exist(isiFile,'file') && exist('ISIspar','var')
    fprintf(1,'Saving the inter-spike intervals for each cluster... ');
    save(isiFile,'ISIspar','ISIVals','-v7.3')
    fprintf(1,'Done!\n')
end
% Number of clusters + the piezo as the first event + the laser as the last
% event, number of time samples in between the time window, and number of
% total triggers.
[Ne, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);
% Select the onset or the offset of a trigger
fprintf(1,'Condition(s):\n')
fprintf('- ''%s''\n', consCondNames)

% Subscript to indicate the conditions with all whisker stimulations, and
% combinations
allWhiskerStimulus = chCond;
Nccond = length(consCondNames);
%% Conditions' boolean flags
delayFlags = false(NTa,Nccond);
counter2 = 1;
for ccond = consCondSubs(:)'
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

% Results directory. Not the best name, but works for now...
resDir = fullfile(dataDir, 'Results');
if ~exist(resDir, 'dir')
    if ~mkdir(resDir)
        fprintf(1, "There was an issue creating %s!\n", resDir)
        fprintf(1, "Saving results in main directory")
        resDir = dataDir;
    end
end
resPttrn = 'Res VW%.2f-%.2f ms %s ms %s ms %s.mat';
resFN = sprintf(resPttrn, timeLapse*1e3, RW_key, SW_key, C_key);
resFP = fullfile(resDir, resFN);

% Statistical scatter figure names
prmSubs = nchoosek(1:Nccond,2); Nsf = size(prmSubs,1) + Nccond;
snglSubs = cumsum(Nccond:-1:1); cmbSubs = setdiff(1:Nsf, snglSubs);
cmpCondNames = string(consCondNames(:));
cmpCondNames = cat(1, cmpCondNames, arrayfun(@(x) ...
    consCondNames(prmSubs(x,1)) + " vs. " + ...
    consCondNames(prmSubs(x,2)), (1:size(prmSubs, 1))'));
stFigSubfix = "";
if metaNameFlag
    stFigSubfix = stFigSubfix + " " + RW_key + " " + SW_key;
end
stFigFN = fullfile(ephFigDir, "Stat " + cmpCondNames + stFigSubfix);
cmpCondNames_aux([snglSubs, cmbSubs]) = stFigFN; stFigFN = cmpCondNames_aux;

if exist(resFP,"file") && all(arrayfun(@(x) exist(x, "file"), stFigFN))
    load(resFP, "Results", "Counts")
    arrayfun(@(x) uiopen(x + ".fig", true), stFigFN)
else
    % Statistical tests
    [Results, Counts] = statTests(discStack, delayFlags, timeFlags);
    % Plotting statistical tests
    [Figs, Results] = scatterSignificance(Results, Counts, consCondNames,...
        delta_t, gclID); configureFigureToPDF(Figs);
    arrayfun(@(x,y) saveFigure(x,y,true), Figs, stFigFN(:))
    save(resFP, "Results", "Counts", "configStructure", "gclID")
end
[rclIdx, H, zH] = getSignificantFlags(Results);
Htc = sum(H,2);
CtrlCond = contains(consCondNames,'control','IgnoreCase',true);
if ~nnz(CtrlCond)
    CtrlCond = true(size(H,2),1);
end
wruIdx = all(H(:,CtrlCond),2);
Nwru = nnz(wruIdx);
fprintf('%d responding clusters:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})

%% Map prototype
mapPttrn = "Map %s.mat";
resMap_path = fullfile(resDir, sprintf(mapPttrn, expName));
nmFlag = true;
try
    resMap = MapNested();
catch
    fprintf(1, "'RolandRitt/Matlab-NestedMap' toolbox not installed! ")
    fprintf(1, "Cannot create multi-key map!\n")
    nmFlag = false;
end
if nmFlag
    if exist(resMap_path, "file")
        resMap_vars = load(resMap_path,"resMap", "keyCell");
        resMap = resMap_vars.resMap;
        keyCell = resMap_vars.keyCell; clearvars respMap_vars
        try
            % Testing if the current key combination exists
            resMap_value = resMap(currentMapKey{:});
            clearvars resMap_value
        catch
            fprintf(1, "Saving responsive unit flags\n")
            resMap(currentMapKey{:}) = wruIdx;
            keyCell = cat(1, keyCell, currentMapKey);
        end
    else
        fprintf(1, "Saving responsive unit flags\n")
        resMap(currentMapKey{:}) = wruIdx;
        keyCell = currentMapKey;
    end
    save(resMap_path, "resMap", "keyCell")
end

%% Filter question
filterIdx = true(Ne,1);
ansFilt = questdlg('Would you like to filter for significance?','Filter',...
    'Yes','No','Yes');
filtStr = 'unfiltered'; filtFlag = false;
if strcmp(ansFilt,'Yes')
    filterIdx = [true; wruIdx]; filtFlag = true;
    filtStr = 'filtered';
end
%% Getting the relative spike times for the whisker responsive units (wru)
% For each condition, the first spike of each wru will be used to compute
% the standard deviation of it.

% cellLogicalIndexing = @(x,idx) x(idx);
isWithinResponsiveWindow =...
    @(x) x > responseWindow(1) & x < responseWindow(2);
relSpkFN = string(expName) + " " + RW_key + " " + SW_key + " " + ...
    VW_key + " ms " + C_key + " (" + string(filtStr) + ") RelSpkTms.mat";
consVars = {'relativeSpkTmsStruct', 'firstSpkStruct', ...
    'SpontaneousStruct', 'configStructure'}; 
rspMF = matfile(fullfile(dataDir, relSpkFN)); 
%{
relSpkFN =...
        sprintf('%s RW%.2f - %.2f ms SW%.2f - %.2f ms VW%.2f - %.2f ms %s (%s) exportSpkTms.mat',...
        expName, responseWindow*1e3, spontaneousWindow*1e3,...
        timeLapse*1e3, Conditions(chCond).name, filtStr);
%}
if ~exist(relSpkFN,'file') || any(~contains(who(rspMF), consVars))
    rst = arrayfun(@(x) getRasterFromStack(discStack, ~delayFlags(:,x), ...
        filterIdx(3:end), timeLapse, fs, true, true), 1:size(delayFlags,2), ...
        fnOpts{:});
    relativeSpkTmsStruct = struct('name', cellstr(consCondNames), ...
        'SpikeTimes', rst);
    firstSpkStruct = getFirstSpikeInfo(relativeSpkTmsStruct, configStructure);
    
    % Spontaneous firing rates
    Texp = Ns/fs;
    trainDuration = 1;
    AllTriggs = unique(cat(1, Conditions.Triggers), 'rows', 'sorted');
    [spFr, ~, SpSpks, spIsi] = getSpontFireFreq(spkSubs, AllTriggs,...
        [0, Texp], fs, trainDuration + delta_t + responseWindow(1));
    SpontaneousStruct = struct('Spikes', SpSpks, 'FR', ...
        arrayfun(@(x) {x}, spFr), 'ISI', spIsi);
    save(fullfile(dataDir, relSpkFN), consVars{:})
else
    load(fullfile(dataDir, relSpkFN), consVars{1:3})
end

%% Ordering PSTH

orderedStr = 'ID ordered';
dans = questdlg('Do you want to order the PSTH other than by IDs?',...
    'Order', 'Yes', 'No', 'No');
ordSubs = 1:nnz(filterIdx(2:Ncl+1));
pclID = gclID(filterIdx(2:Ncl+1));
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
        [~,ordSubs] = sortrows(clInfo(pclID,:),ordVar);
    end
end

%% Plot PSTH
if exist('Triggers', 'var')
    csNames = fieldnames(Triggers);
end
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end

%psthTx = (0:Nbn-1) * binSz + timeLapse(1);
 Ntc = size(cst,2);
psthFN = "PSTH " + consCondNames(:) + " " + BZ_key + " " + string(orderedStr);
if filtFlag
    psthFN = psthFN + " " + filtStr;
end
% PSTH construction
psthFP = fullfile(ephFigDir, psthFN);
[PSTH, trig] = arrayfun(@(x) getPSTH(discStack(filterIdx,:,:), ...
    timeLapse, ~delayFlags(:,x), binSz, fs), 1:Nccond, fnOpts{:});
if any(arrayfun(@(x) ~exist(x+".fig","file"), psthFP))
    if exist('cst', 'var') && ~isempty(cst)
        % Take into account covariance for signals.
        stims = arrayfun(@(x) mean(cst(:,:,delayFlags(:,x)),3), 1:Nccond, ...
            fnOpts{:}); stims = cellfun(@(x) zscore(x, 0, 'all'), stims, fnOpts{:});
    else
        stims = repmat({zeros(1,Ntc)}, Nccond, 1);
    end
    psthFigs = cellfun(@(p,t,n,ids,s) plotClusterReactivity(p(ordSubs,:), t,...
        n, timeLapse, binSz, [ids; pclID(ordSubs)], strrep(expName,'_',' '), ...
        s, csNames), PSTH, trig, num2cell(Na), cellstr(consCondNames), stims);
    arrayfun(@(f) ylabel(f.Children(end), ...
        [f.Children(end).YLabel.String, sprintf('^{%s}',orderedStr)]), ...
        psthFigs);
    arrayfun(@(f, fn) saveFigure(f, fn), psthFigs(:), psthFP(:));
else
   arrayfun(@(f) uiopen(f + ".fig", true), psthFP) 
end
% Z-score PSTH for all units
ephysPttrn = 'Z-score all-units PSTH %s Ntrials%s';
ephysName = sprintf(ephysPttrn, sprintf('%s ', consCondNames{:}), ...
    sprintf(' %d', Na));
ephysFile = fullfile(ephFigDir, ephysName);
if ~exist(ephysFile, 'file')
    [ppFig, PSTHall] = compareCondPSTHs(cat(3,PSTH{:}), Na, binSz, ...
        timeLapse, consCondNames);
    saveFigure(ppFig, ephysFile, 1);
else
    uiopen(ephysFile, true);
end
clearvars ppFig ephysPttrn ephysName ephysFile
%% Log PSTH
Nbin = 64;
ncl = size(relativeSpkTmsStruct(1).SpikeTimes,1);

logPSTH = getLogTimePSTH(relativeSpkTmsStruct, true(ncl,1),...
    'tmWin', responseWindow, 'Offset', 2.5e-3, 'Nbin', Nbin,...
    'normalization', 'fr');
lpFN = sprintf("Log-likePSTH %s %d-conditions NB%d",...
    logPSTH.Normalization, Nccond, Nbin);
if Nccond > 1
    lmiFN = sprintf("LogMI %d-conditions NB%d", Nccond, Nbin);
    lmiFP = fullfile(ephFigDir, lmiFN);
end
if filtFlag
    lpFN = lpFN + " (" + filtStr + ")";
    if Nccond > 1 
        lmiFP = lmiFP + " (" + filtStr + ")";
    end
end
lpFP = fullfile(ephFigDir, lpFN);
if ~exist(lpFP+".fig", "file")
    logFigs = plotLogPSTH(logPSTH); saveFigure(logFigs(1), lpFP, true)
    if numel(logFigs) > 1
        saveFigure(logFigs(2), lmiFP, true)
        popEffects = logFigs(2).UserData; vrs = who(matfile(resFP));
        MIStruct = struct('ConditionNames', consCondNames, ...
            'MI', arrayfun(@(x) struct('Comparative', ...
            string(consCondNames(popEffects(x,1)))+" vs "+...
            string(consCondNames(popEffects(x,2))), 'Value', popEffects(x,3)), ...
            1:size(popEffects,1), fnOpts{:}));
        if ~any(ismember(vrs,'MIStruct'))
            fprintf(1,'Adding "MIStruct" to %s\n', resFN)
            save(resFP, 'MIStruct','-append')
        end
    end
else
    uiopen(lpFP+".fig", true)
    load(resFP, "MIstruct")
    if Nccond > 1
        uiopen(lmiFP+".fig",true)
    end
end
%% Cluster population proportions
% Responsive and unresponsive cells, significantly potentiated or depressed
% and unmodulated.
fnOpts = {'UniformOutput', false};
hsOpts = {'BinLimits', [-1,1], 'NumBins', 32,...
    'Normalization', 'probability', 'EdgeColor', 'none', 'DisplayName'};
axOpts = {'Box', 'off', 'Color', 'none'};
evFr = cellfun(@(x) mean(x,2)./diff(responseWindow), Counts(:,2),fnOpts{:});
evFr = cat(2, evFr{:});
getMI = @(x) diff(x, 1, 2)./sum(x, 2);
MIevok = getMI(evFr);
Nrn = sum(wruIdx); Ntn = size(wruIdx,1);
signMod = Results(1).Activity(2).Pvalues < 0.05;
potFlag = MIevok > 0;
Nrsn = sum(wruIdx & signMod); Nrsp = sum(wruIdx & signMod & potFlag);
% All spikes in a cell format
if size(spkSubs,1) < size(gclID,1)
    spkSubs = cat(1, {round(sortedData{goods(1),2}*fs)}, spkSubs);
end
NaCount = 1;
% Experiment firing rate and ISI per considered condition
spFrC = zeros(Ncl, size(consCondSubs(:),1), 'single');
econdIsi = cell(Ncl, size(consCondSubs(:),1));
econdSpks = econdIsi;
trainDuration = 1;
for ccond = consCondSubs(:)'
    itiSub = mean(diff(Conditions(ccond).Triggers(:,1)));
    consTime = [Conditions(ccond).Triggers(1,1) - round(itiSub/2),...
        Conditions(ccond).Triggers(Na(NaCount),2) + round(itiSub/2)]...
        ./fs;
    [spFrC(:,NaCount),~, econdSpks(:,NaCount), econdIsi(:,NaCount)] =...
        getSpontFireFreq(spkSubs, Conditions(ccond).Triggers,...
        consTime, fs, trainDuration + delta_t + responseWindow(1));
    NaCount = NaCount + 1;
end
MIspon = getMI(spFrC); SNr = evFr./spFrC;
%% Plot proportional pies
clrMap = lines(2); clrMap([3,4],:) = [0.65;0.8].*ones(2,3);
% Responsive and non responsive clusters
respFig = figure("Color", "w");
pie([Ntn-Nrn, Nrn], [0, 1], {'Unresponsive', 'Responsive'});
pObj = findobj(respFig, "Type", "Patch");
arrayfun(@(x) set(x, "EdgeColor", "none"), pObj);
arrayfun(@(x) set(pObj(x), "FaceColor", clrMap(x+2,:)), 1:length(pObj))
propPieFileName = fullfile(ephFigDir,...
    sprintf("Whisker responsive proportion pie %s (%dC, %dR)",...
    C_key, [Ntn-Nrn, Nrn]));
saveFigure(respFig, propPieFileName, 1);
% Potentiated, depressed and unmodulated clusters pie
if Nccond == 2
    potFig = figure("Color", "w");
    pie([Nrn - Nrsn, Nrsp, Nrsn - Nrsp], [0, 1, 1], {'Non-modulated', ...
        'Potentiated', 'Depressed'}); % set(potFig, axOpts{:})
    pObj = findobj(potFig, "Type", "Patch");
    arrayfun(@(x) set(x, "EdgeColor", "none"), pObj);
    arrayfun(@(x) set(pObj(x), "FaceColor", clrMap(x,:)), 1:length(pObj))
    modPropPieFigFileName = fullfile(ephFigDir,...
        sprintf("Modulation proportions pie %s (%dR, %dP, %dD)",...
        C_key, Nrn - Nrsn, Nrsp, Nrsn - Nrsp));
    saveFigure(potFig, modPropPieFigFileName, 1)
    % Modulation index histogram
    MIFig = figure; histogram(MIspon, hsOpts{:}, "Spontaneous"); hold on;
    histogram(MIevok, hsOpts{:}, "Evoked"); set(gca, axOpts{:});
    title("Modulation index distribution"); xlabel("MI");
    ylabel("Cluster proportion"); lgnd = legend("show");
    set(lgnd, "Box", "off", "Location", "best")
    saveFigure(MIFig, fullfile(ephFigDir,...
        "Modulation index dist evoked & after induction "+C_key), 1)
end
%% Get significantly different clusters
gcans = questdlg(['Do you want to get the waveforms from the',...
    ' ''responding'' clusters?'], 'Waveforms', 'Resp only', 'All', 'None', ...
    'None');
switch gcans
    case 'Resp only'
        clWaveforms = getClusterWaveform(gclID(wruIdx), dataDir);
    case 'All'
        clWaveforms = getClusterWaveform(gclID, dataDir);
    case 'None'
        fprintf(1, 'You can always get the waveforms later\n')
end

%% Rasters from interesting clusters
rasAns = questdlg('Plot rasters?','Raster plot','Yes','No','Yes');
if strcmpi(rasAns,'Yes')
    [clSel, iOk] = listdlg('ListString',pclID,...
        'Name', 'Selection of clusters',...
        'CancelString', 'None',...
        'PromptString', 'Select clusters',...
        'SelectionMode', 'multiple');
    if ~iOk
        return
    end
    % Color of the rectangle
    cmap = [{[0.6, 0.6, 0.6]};... gray
        {[1, 0.6, 0]};... orange
        {[0.5, 0.6, 0.5]};... greenish gray (poo)
        {[0, 0.6, 1]};... baby blue
        {[0.6, 0.6, 1]};... lila
        {[1, 0.6, 1]}]; % pink
    clrNames = {'Gray','Orange','Poo','Baby Blue','Lila','Pink'}';
    clrmap = containers.Map(clrNames,cmap);
    if contains(Conditions(chCond).name, whStim,'IgnoreCase', 1)
        rectColor = clrmap('Orange');
    elseif contains(Conditions(chCond).name, cxStim, 'IgnoreCase', 1)
        rectColor = clrmap('Baby Blue');
    else
        rectColor = clrmap('Gray');
    end
    % Choose the conditions to be plotted
    resCondNames = arrayfun(@(x) x.name, Conditions(consCondSubs),...
        'UniformOutput', 0);
    [rasCondSel, iOk] = listdlg('ListString', resCondNames,...
        'PromptString', 'Which conditions to plot?',...
        'InitialValue', 1:length(consCondSubs),...
        'CancelString', 'None',...
        'Name', 'Conditions for raster',...
        'SelectionMode', 'multiple');
    if ~iOk
        return
    end
    rasCond = consCondSubs(rasCondSel);
    rasCondNames = consCondNames(rasCondSel);
    Nrcl = numel(clSel);
    % Reorganize the rasters in the required order.
    clSub = find(ismember(gclID, pclID(clSel)))+1;
    [rasIdx, rasOrd] = ismember(pclID(ordSubs), pclID(clSel));
    clSub = clSub(rasOrd(rasIdx));
    clSel = clSel(rasOrd(rasOrd ~= 0));
    Nma = min(Na(rasCondSel));
    rasFig = figure;
    Nrcond = length(rasCond);
    ax = gobjects(Nrcond*Nrcl,1);
    timeFlags = all([tx(:) >= timeLapse(1), tx(:) <= timeLapse(2)],2);
    for cc = 1:length(rasCond)
        % Equalize trial number
        trigSubset = sort(randsample(Na(rasCondSel(cc)),Nma));
        tLoc = find(delayFlags(:,rasCondSel(cc)));
        tSubs = tLoc(trigSubset);
        % Trigger subset for stimulation shading
        trigAlSubs = Conditions(rasCond(cc)).Triggers(trigSubset,:);
        timeDur = round(diff(trigAlSubs, 1, 2)/fs, 3); tol = 0.02/max(timeDur);
        timeDurUniq = uniquetol(timeDur, tol);
        %trigChange = find(diff(timeDur) ~= 0);
        trigChange = find(~ismembertol(timeDur, timeDurUniq, tol));
        for ccl = 1:Nrcl
            lidx = ccl + (cc - 1) * Nrcl;
            ax(lidx) = subplot(Nrcond, Nrcl, lidx);
            title(ax(lidx),sprintf('%s cl:%s',rasCondNames{cc},pclID{clSel(ccl)}))
            plotRasterFromStack(discStack([1,clSub(ccl)],timeFlags,tSubs),...
                timeLapse, fs,'',ax(lidx));
            ax(lidx).YAxisLocation = 'origin';ax(lidx).YAxis.TickValues = Nma;
            ax(lidx).YAxis.Label.String = Nma;
            ax(lidx).YAxis.Label.Position =...
                [timeLapse(1)-timeLapse(1)*0.65, Nma,0];
            %ax(lidx).XAxis.TickLabels =...
            %    cellfun(@(x) str2num(x)*1e3, ax(lidx).XAxis.TickLabels,...
            %    'UniformOutput', 0);
            xlabel(ax(lidx), 'Time [s]')
            initSub = 0;
            optsRect = {'EdgeColor',rectColor,'FaceColor','none'};
            for ctr = 1:numel(trigChange)
                rectangle('Position',[0, initSub,...
                    timeDur(trigChange(ctr)), trigChange(ctr)],optsRect{:})
                initSub = trigChange(ctr);
            end
            rectangle('Position', [0, initSub, timeDur(Nma),...
                Nma - initSub],optsRect{:})
        end
    end
    linkaxes(ax,'x')
    rasFigName = sprintf('%s R-%scl_%sVW%.1f-%.1f ms', expName,...
        sprintf('%s ', rasCondNames{:}), sprintf('%s ', pclID{clSel}),...
        timeLapse*1e3);
    rasFigPath = fullfile(ephFigDir, rasFigName);
    arrayfun(@(x) set(x,'Color','none'), ax);
    saveFigure(rasFig, rasFigPath, 1);
    clearvars ax rasFig
end
%% Response speed characterization
btx = (0:Nbn-1)*binSz + timeLapse(1);
respIdx = isWithinResponsiveWindow(btx);
% Window defined by the response in the population PSTH
% twIdx = btx >= 2e-3 & btx <= 30e-3;
% respTmWin = [2, 30]*1e-3;
[mdls, r2, qVals, qDiff] = exponentialSpread(PSTH{1}, btx, responseWindow);
mdls(mdls(:,2) == 0, 2) = 1;
%% Cross-correlations
ccrAns = questdlg(['Get cross-correlograms?',...
    '(Might take a while to compute if no file exists!)'],...
    'Cross-correlations', 'Yes', 'No', 'Yes');
if strcmpi(ccrAns,'Yes')
    consTime = inputdlg('Cross-correlation lag (ms):',...
        'Correlation time', [1,30], {'50'});
    try
        consTime = str2num(consTime{1})*1e-3;
    catch
        fprintf(1, 'Cancelling...\n');
        return
    end
    corrsFile = fullfile(dataDir,sprintf('%s_%.2f ms_ccorr.mat', expName,...
        consTime*1e3));
    if ~exist(corrsFile,'file')
        if size(spkSubs,1) < size(goods,1)
            spkSubs = cat(1, {round(sortedData{goods(1),2}*fs)},spkSubs);
        end
        corrs = neuroCorr(spkSubs, consTime, 1, fs);
        save(corrsFile,'corrs','consTime','fs','-v7.3');
    else
        load(corrsFile, 'corrs')
    end
end
%% Auto-correlation measurement
% Arranging the auto-correlograms out of the cross-correlograms into a
% single matrix
try
    acorrs = cellfun(@(x) x(1,:), corrs, 'UniformOutput', 0);
catch
    fprintf(1, 'No correlograms in the workspace!\n')
end
if exist('acorrs', "var")
    acorrs = single(cat(1, acorrs{:})); Ncrs = size(acorrs, 2);
    % Computing the lag axis for the auto-correlograms
    b = -ceil(Ncrs/2)/fs; corrTx = (1:Ncrs)'/fs + b;
    corrTmWin = [1, 25]*1e-3;
    % Interesting region in the auto-correlogram
    [cmdls, cr2, cqVals, cqDiff] =...
        exponentialSpread(acorrs, corrTx, corrTmWin);
end
% lwIdx = corrTx >= (1e-3 - 1/fs) & corrTx <= 25e-3;
% lTx = corrTx(lwIdx);
% icsAc = double(1 - cumsum(acorrs(:,lwIdx)./sum(acorrs(:,lwIdx),2),2));
% % Quartile cuts
% quartCut = exp(-log([4/3, 2, exp(1), 4]))';
% % Exponential analysis for the auto-correlograms
% cmdls = zeros(Ncl,2); cr2 = zeros(Ncl,1); cqVals = zeros(Ncl,4);
% for ccl = 1:Ncl
%     % Exponential fit for the inverted cumsum
%     [fitObj, gof] = fit(lTx, icsAc(ccl,:)', 'exp1','Display','off');
%     cmdls(ccl,:) = coeffvalues(fitObj); cr2(ccl) = gof.rsquare;
%     % Quartiles cut for exponential distribution (25, 50, 63.21, 75)
%     quartFlags = icsAc(ccl,:) >= quartCut;
%     [qSubs, ~] = find(diff(quartFlags'));
%     il = arrayfun(@(x) fit_poly(lTx(x:x+1), icsAc(ccl,x:x+1), 1), qSubs,...
%         'UniformOutput', 0);il = cat(2,il{:});
%     cqVals(ccl,:) = (quartCut' - il(2,:))./il(1,:);
% end
% cqDiff = diff(cqVals(:,[1,4]),1,2);

%% Behaviour

flds = dir(getParentDir(dataDir,1));
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), flds);
flds(pointFlag) = [];
behFoldFlag = arrayfun(@(x) any(strcmpi(x.name, 'Behaviour')), flds);
possNames = ["P", "L"];
if any(behFoldFlag) && sum(behFoldFlag) == 1
    % If only one folder named Behaviour exists, chances are that this is
    % an awake experiment.
    behDir = fullfile(flds(behFoldFlag).folder,flds(behFoldFlag).name);
    fprintf(1, "Found %s!\n", behDir)
    answ = questdlg('Analyse behaviour?','Behaviour','Yes','No','Yes');
    if strcmpi(answ,'Yes')
        hstOpts = {'BinMethod', 'integers', 'BinLimits', [-0.5,4.5]};
        behChCond = cellfun(@(x) contains(Conditions(chCond).name, x), ...
            {["Piezo", "Puff"];["Laser","Light"]});
        [behRes, behFigDir] = analyseBehaviour(behDir, 'Condition', possNames(behChCond), ...
            'PairedFlags', delayFlags, 'FigureDirectory', FigureDir, ...
            'ConditionsNames', cellstr(consCondNames));
        biFigPttrn = "BehIndex%s";
        biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames));
        [pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
        behRes = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), behRes, pAreas);
        set(behAreaFig, 'UserData', behRes)
        
        biFN = sprintf(biFigPttrn, pAreas);
        saveFigure(behAreaFig, fullfile(behFigDir, biFN), true, true);

        trMvFlag = arrayfun(@(cr) behRes(1).Results(cr).MovStrucure.MovmentFlags, ...
            1:size(behRes(1).Results,2), fnOpts{:}); trMvFlag = cat(3, trMvFlag{:});
        BIscaleMat = sum(trMvFlag,3);
        BIscale = arrayfun(@(cc) BIscaleMat(pairedStim(:,cc), cc), 1:Nccond, ...
            fnOpts{:});
        [hg, hg_bin] = cellfun(@(c) histcounts(c, hstOpts{:}), ...
            BIscale, fnOpts{:}); 
        hg = cat(1, hg{:}); hg_bin = cat(1, hg_bin{:});
        %%
        figure; bar(0:4, hg./sum(hg,2), 'EdgeColor', 'none'); hold on;
        ylim([0,1]); set(gca, axOpts{:}); 
        legend(consCondNames, 'AutoUpdate','off', lgOpts{:})
        lmbdaHeight = 0.95-(0.15/Nccond)*(0:Nccond-1);
        arrayfun(@(pd) scatter(poaDist(pd).lambda, lmbdaHeight(pd), '|',...
            'MarkerEdgeColor', clrMap(pd,:)), 1:Nccond)
        
        arrayfun(@(pd) line(paramci(poaDist(pd)), ...
            lmbdaHeight([pd,pd]), 'Color', clrMap(pd,:), ...
            'Marker', '|'), 1:Nccond)
        title(strrep(expName, '_',' ')); xlabel('Moving body parts')
        ylabel('Trial proportion')
    end
end