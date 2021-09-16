% 30.08.19 Jittering analysis by using the Data Explorer. 
clearvars
%% Load the data
% Choosing the working directory
dataDir = uigetdir('E:\Data\VPM\Jittering\Silicon Probes\',...
    'Choose a working directory');
if dataDir == 0
    return
end
% Creating the figure directory
figureDir = fullfile(dataDir,'Figures\');
if ~mkdir(figureDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end
% Loading the necessary files
if ~loadTriggerData(dataDir)
    fprintf(1,'Not possible to load all the necessary variables\n')
    return
end

%% Constructing the helper 'global' variables
% Number of total samples
Ns = structfun(@numel,Triggers);
Ns = min(Ns(Ns>1));
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
    clInfo = addvars(clInfo,~badsIdx,'After','id',...
        'NewVariableNames','ActiveUnit');
    writeClusterInfo(clInfo, fullfile(dataDir, 'cluster_info.tsv'), 1);
end
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
%% Inter-spike intervals
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
    try
        ISIspar = sparse(rows, cols, vals);
    catch
        fprintf(1, 'Not possible to create such a big array')
    end
else
    load(isiFile,'ISIspar')
end
% ISIsignal = zeros(Ncl,Ns,'single');
% for ccl = 1:Ncl
%     ISIsignal(ccl,spkSubs2{ccl}) = ISIVals{ccl};
% end
%% User controlling variables
% Time lapse, bin size, and spontaneous and response windows
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
fprintf(1,'Time window: %.2f - %.2f ms\n',timeLapse(1)*1e3, timeLapse(2)*1e3)
fprintf(1,'Response window: %.2f - %.2f ms\n',responseWindow(1)*1e3, responseWindow(2)*1e3)
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

%% Constructing the stack out of the user's choice
% discStack - dicrete stack has a logical nature
% cst - continuous stack has a numerical nature
% Both of these stacks have the same number of time samples and trigger
% points. They differ only in the number of considered events.
[discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
    timeLapse,fs,fs,spkSubs,continuousSignals);
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
    saveFigure(Figs(cc), stFigName)
end
H = cell2mat(cellfun(@(x) x.Pvalues,...
    arrayfun(@(x) x.Activity, Results(indCondSubs), 'UniformOutput', 0),...
    'UniformOutput', 0)) < 0.05;

Htc = sum(H,2);
CtrlCond = contains(consCondNames,'control','IgnoreCase',true);
if ~nnz(CtrlCond)
    CtrlCond = true(size(H,2),1);
end
wruIdx = all(H(:,CtrlCond),2);
Nwru = nnz(wruIdx);

fprintf('%d whisker responding clusters:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})
%% Get significantly different clusters
gcans = questdlg(['Do you want to get the waveforms from the',...
    ' ''responding'' clusters?'], 'Waveforms', 'Yes', 'No', 'No');
if strcmp(gcans, 'Yes')
    clWaveforms = getClusterWaveform(gclID(wruIdx), dataDir);
end

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
            'Mechanical',reshape(mean(cst(whFlag,:,delayFlags(:,cdel)),3),...
            1,Nt),'Laser',reshape(mean(cst(lrFlag,:,delayFlags(:,cdel)),3),...
            1,Nt),'TimeAxis',(0:Nt-1)/fs + timeLapse(1));
        cdel = cdel + 1;
    end
    save(fullfile(dataDir,[expName,'analysis.mat']),'Conditions','-append')
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
    filtStr = 'filtered';
end

%% Getting the relative spike times for the whisker responsive units (wru)
% For each condition, the first spike of each wru will be used to compute
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
csvSubfx = sprintf(' VW%.1f-%.1f ms (%s).csv', timeLapse*1e3, filtStr);
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
%     respIdx = cellfun(isWithinResponsiveWindow, relativeSpikeTimes,...
%         'UniformOutput',false);
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
    if ~existFlag
        fID = fopen(csvFileName,'w');
        fprintf(fID,'%s, %s\n','Cluster ID','Relative spike times [ms]');
    end
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
    %{ 
    % First spike in the trial
    spikeTimesINRespWin = cellfun(cellLogicalIndexing,...
        relativeSpikeTimes, respIdx, 'UniformOutput',false);
    allSpikeTimes = cell2mat(spikeTimesINRespWin(:)');
    condParams(:,:,ccond) = emforgmm(allSpikeTimes, M, 1e-6, 0);
    condPDF(:,ccond) = genP_x(condParams(:,:,ccond), txpdf);
    firstOrdStats(:,ccond) = [mean(allSpikeTimes), std(allSpikeTimes)];
    hfig = figure('Visible', 'off'); h = histogram(allSpikeTimes, binAx,...
        'Normalization', 'probability');
    condHist(:,ccond) = h.Values;
    close(hfig)
    for ccl = 1:Nwru
        frstSpikeFlag = ~cellfun(@isempty,spikeTimesINRespWin(ccl,:));
        firstSpike(ccl,ccond) = std(...
            cell2mat(spikeTimesINRespWin(ccl,frstSpikeFlag)));    
    end
    %}
end

relSpkFileName =...
    sprintf('%s RW%.2f - %.2f ms SW%.2f - %.2f ms %s (%s) exportSpkTms.mat',...
    expName, responseWindow*1e3, spontaneousWindow*1e3,...
    Conditions(chCond).name, filtStr);
if ~exist(relSpkFileName,'file')
    save(fullfile(dataDir, relSpkFileName), 'relativeSpkTmsStruct',...
        'configStructure')
end
%% Standard Deviations of First Spikes After Each Trigger per Unit
% firstSpikes(relativeSpkTmsStruct, gclID, dataDir);
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
goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
psthTx = (0:Nbn-1) * binSz + timeLapse(1);
psthFigs = gobjects(Nccond,1);
for ccond = 1:Nccond
    figFileName =...
        sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, consCondNames{ccond}, timeLapse*1e3, binSz*1e3,...
        responseWindow*1e3, spontaneousWindow*1e3, onOffStr, orderedStr,...
        filtStr);
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
    psthFigs(ccond) = plotClusterReactivity(PSTH(ordSubs,:,ccond), trig,...
        sweeps, timeLapse, binSz, [consCondNames(ccond); pclID(ordSubs)],...
        strrep(expName,'_','\_'), stims, csNames);
    psthFigs(ccond).Children(end).YLabel.String =...
        [psthFigs(ccond).Children(end).YLabel.String,...
        sprintf('^{%s}',orderedStr)];
    figFilePath = fullfile(figureDir, figFileName);
    saveFigure(psthFigs(ccond), figFilePath);
end
%% Log PSTH -- Generalise this part!!
Nbin = 64;
ncl = size(relativeSpkTmsStruct(1).SpikeTimes,1);


logPSTH = getLogTimePSTH(relativeSpkTmsStruct, true(ncl,1),...
    'tmWin', responseWindow, 'Offset', 2.5e-3, 'Nbin', Nbin,...
    'normalization', 'fr');
logFigs = plotLogPSTH(logPSTH);
% Saving the figures
lpFigName = sprintf('%s Log-likePSTH %s %d-conditions RW%.1f-%.1f ms NB%d (%s)',...
    expName, logPSTH.Normalization, Nccond, responseWindow*1e3, Nbin, filtStr);
saveFigure(logFigs(1), fullfile(figureDir, lpFigName))
if numel(logFigs) > 1
    lmiFigName = sprintf('%s LogMI %d-conditions RW%.1f-%.1f ms NB%d (%s)',...
        expName, Nccond, responseWindow*1e3, Nbin, filtStr);
    saveFigure(logFigs(2), fullfile(figureDir, lmiFigName))
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
    resCondNames = arrayfun(@(x) x.name, Conditions(consideredConditions),...
        'UniformOutput', 0);
    [rasCondSel, iOk] = listdlg('ListString', resCondNames,...
        'PromptString', 'Which conditions to plot?',...
        'InitialValue', 1:length(consideredConditions),...
        'CancelString', 'None',...
        'Name', 'Conditions for raster',...
        'SelectionMode', 'multiple');
    if ~iOk
        return
    end
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
    rasFigPath = fullfile(figureDir, rasFigName);
    arrayfun(@(x) set(x,'Color','none'), ax); 
    saveFigure(rasFig, rasFigPath);
    clearvars ax rasFig
end
%% Response speed characterization
btx = (0:Nbn-1)*binSz + timeLapse(1);
respIdx = isWithinResponsiveWindow(btx);
% Window defined by the response in the population PSTH
% twIdx = btx >= 2e-3 & btx <= 30e-3;
% respTmWin = [2, 30]*1e-3;
[mdls, r2, qVals, qDiff] = exponentialSpread(PSTH(:,:,1), btx, responseWindow);
mdls(mdls(:,2) == 0, 2) = 1;
%% Optotag
% Clusters with 0.63 spikes per stimulus in a 1 ms bin will be considered
% with consistency, 50% of the response before 7 ms as readily available,
% and with a small spread less than 2 ms in between 1 and 3 quartile as
% precise. Considering unfiltered PSTH and the spikes within the response
% window. Experiments with ChR2
rangeAll = @(x) [min(x), max(x)];
if contains(Conditions(chCond).name,'laser','IgnoreCase',1)
    fprintf(1,'Optotagging clusters... ')
    tvNames = clInfo.Properties.VariableNames;
    optoVarFlag = contains(tvNames, 'Optotag');
    prevFlag = false;
    if any(optoVarFlag)
        fprintf(1, 'Previous optotag found!\n')
        optoIdx = clInfo{clInfo.ActiveUnit == 1, 'Optotag'} == 1;
        prevFlag = true;
    else
        PSTHtrial = PSTH ./ Na; 
        
        
        oqVals = qVals(optoCl,:);
        [~, modeTm] = max(optoPSTH(:, respIdx, :),[],2);
        availableIdx = oqVals(:,3) <= 7e-3; % Availability
        preciseIdx = (oqVals(:,5) - oqVals(:,2)) < 2e-3; % Precision
        optoTaggedCl = optoCl(availableIdx & preciseIdx);
        optoIdx = ismember(gclID,pclID(optoTaggedCl));
        fprintf(1, 'Found %d clusters\n', sum(optoIdx))
    end
    if any(optoIdx)
        
        % Probe view
        chPos = readNPY(fullfile(dataDir, 'channel_positions.npy'));
        chMap = readNPY(fullfile(dataDir, 'channel_map.npy')); scl = 10;
        mdX = mean(rangeAll(chPos(:,1)));
        % Matrix for scaling the xlimit for the probe view.
        lsrCMap = [0.2, 0.8, 1];
        optoOpts = {'MarkerEdgeColor',lsrCMap,'MarkerEdgeAlpha',0.3};
        optoLineOpts = {'Marker','+', 'DisplayName','Optotagged region',...
            'Color',lsrCMap}; scaleMatrix = eye(2)+[-scl;scl].*flip(eye(2),1);
        xview = scaleMatrix * repmat(mdX,2,1);
        probFig = figure('Name', 'Probe', 'Color', [1,1,1]);
        probAx = axes('NextPlot', 'add'); xlim(probAx, xview');
        probPts = scatter(probAx, chPos(:,1), chPos(:,2), '.k');
        title(probAx, 'Electrodes position in the probe');
        probAx.XAxis.Visible = 'off'; text(mdX, -1, 'Tip', 'Parent', probAx,...
            'HorizontalAlignment', "center", "VerticalAlignment",  "top")
        tvNames = clInfo.Properties.VariableNames; tvNames = string(tvNames);
        chanStr = ["ch";"channel"]; [~, varSel] = find(tvNames == chanStr);
        optoCh = clInfo{gclID(optoIdx), varSel}';
        [coordSubs, ~] = find(chMap == optoCh);
        optoPts = scatter(probAx, chPos(coordSubs,1), chPos(coordSubs,2),...
            optoOpts{:}); optoRange = rangeAll(chPos(coordSubs,2));
        optoLine = line(probAx, [mdX;mdX], optoRange', optoLineOpts{:});
        ylabel(probAx, 'Electrode position relative to the probe tip [mm]')
        fprintf(1, 'Optotagged clusters found between %.1f and %.1f mm',...
            optoRange); fprintf(1, ' relative to the tip of the probe.\n')
        legend(probAx, optoLine); saveFigure(probFig, fullfile(figureDir,...
            'Probe view + position of optotagged clusters'))
        if ~prevFlag
            clInfo = addvars(clInfo, false(size(clInfo,1),1),...
                'NewVariableNames', 'Optotag');
            clInfo{gclID(optoIdx), 'Optotag'} = true;
            writeClusterInfo(clInfo, fullfile(dataDir, 'cluster_info.tsv'),...
                1)
        end
    else
        fprintf(1, 'No optotagged clusters found!\n')
    end
end
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
        spkSubs = cat(1, {round(sortedData{goods(1),2}*fs)},spkSubs);
        corrs = neuroCorr(spkSubs, consTime, 1, fs);
        save(corrsFile,'corrs','consTime','fs');
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
    return
end
acorrs = single(cat(1, acorrs{:})); Ncrs = size(acorrs, 2);
% Computing the lag axis for the auto-correlograms
b = -ceil(Ncrs/2)/fs; corrTx = (1:Ncrs)'/fs + b;
corrTmWin = [1, 25]*1e-3;
% Interesting region in the auto-correlogram
[cmdls, cr2, cqVals, cqDiff] =...
    exponentialSpread(acorrs, corrTx, corrTmWin);

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


