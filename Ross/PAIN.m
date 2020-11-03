% PAIN
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
badsIdx = badsIdx | silentUnits;
if ~any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
    clInfo = addvars(clInfo,~badsIdx,'After','id',...
        'NewVariableNames','ActiveUnit');
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
    ISIspar = sparse(rows, cols, vals);
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
defInputs = {'-2, 6', '0.1, 1.5', '0.01'};
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
    fprintf(1, 'Perhaps there''s not enough memory to cope with this ISIs\n');
    % General ISI
    
end
%[dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
%    'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
%    continuousSignals);
if ~exist(isiFile,'file')
    fprintf(1,'Saving the inter-spike intervals for each cluster... ');
    save(isiFile,'ISIspar','isiStack','-v7.3')
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

fprintf('%d whisker responding clusters:\n', Nwru);
fprintf('- %s\n',gclID{wruIdx})
%% Getting cluster info and adding variables to table
% clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
% chanMap = readNPY(fullfile(dataDir,'channel_map.npy'));
% chanPos = readNPY(fullfile(dataDir,'channel_positions.npy'));
% [m,b] = lineariz(chanPos(:,1), 6, 1);
% shank = m*chanPos(:,1) + b;
% shank = round(shank);
% shMap = containers.Map(chanMap, shank);
% setShank = @(x) shMap(x);
% clInfo.shank = arrayfun(setShank, clInfo.channel);
% tb = size(clInfo);
% sz = tb(1);
% ActiveUnit = false(sz,1);
% clInfo = addvars(clInfo,ActiveUnit,'NewVariableNames','ActiveUnit','After','id');
% clInfo{gclID, 'ActiveUnit'} = true;

for a = 1: length(consCondNames)
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_Counts_Spont']} = mean(Counts{a,1}')';
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_Counts_Evoked']} = mean(Counts{a,2}')';
end


% Significant mechanical responses per condition
b = length(consCondNames);
for a = 1 : length(consCondNames)
    clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a}, '_MR']} = Results(b).Activity(1).Pvalues < 0.05;
    b = b + length(consCondNames) - a;
end

% Comparing Across Conditions (i.e. manipulation effect on spontaneous and evoked activity)
% RAM: comment to explain what you are doing
sZ = size(Counts);
d = length(consCondNames);
e = 1;
f = length(consCondNames);
for a = 1:length(consCondNames) - 1
    for b = (a + 1): length(consCondNames)
        
        for c = 1: sZ(1,2)
            clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a},'_vs_',consCondNames{1,b}, '_', Results(e).Activity(c).Type, '_Response']} = Results(e).Activity(c).Pvalues < 0.05;
            %RAM example
            %myString = [consCondNames{1,a},'_vs_',consCondNames{1,b}, '_', Results(e).Activity(c).Type, '_Response'];
            %clInfo{clInfo.ActiveUnit == true,myString} = Results(e).Activity(c).Pvalues < 0.05;
            
        end
        
        e = e + 1;
        
        if e == f
            e = e + 1;
        end
    end
    d = d - 1;
    f = f + d;
end
writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));
%% Get significantly different clusters
gcans = questdlg(['Do you want to get the waveforms ofr the good clusters?'], 'Waveforms', 'Yes', 'No', 'No');
if strcmp(gcans, 'Yes')
    clWaveforms = getClusterWaveform(gclID, dataDir);
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
save(fullfile(dataDir,[expName,'_exportSpkTms.mat']),...
    'relativeSpkTmsStruct','configStructure')
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
    if ~exist([figFileName,'.pdf'], 'file')
        print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.pdf']),...
            '-dpdf','-fillpage')
    end
    if ~exist([figFileName,'.emf'], 'file')
        print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.emf']),...
            '-dmeta')
    end
    
end
for a = 1:length(consideredConditions)
    savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH.fig']));
end
%% Rasters
% DE_Jittering needs to be unfiltered for significance for this to work!


rasterDir = fullfile(figureDir,'Rasters\');
if ~mkdir(rasterDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end
for pwr = [5, 10, 15]
MchTblInd = ['Mech_Control_', num2str(pwr), 'mW_MR'];
LasTblInd = ['Laser_Control_', num2str(pwr), 'mW_LR'];
MchCondControl = ['Mech_Control_', num2str(pwr), 'mW'];
LasCondControl = ['Laser_Control_', num2str(pwr), 'mW'];
MchLasCond = ['Mech_Laser_', num2str(pwr), 'mW'];
EffectTblInd = ['Mech_Control_', num2str(pwr), 'mW_vs_Mech_Laser_', num2str(pwr), 'mW_Evoked_Response'];
TblInd = find(clInfo.(MchTblInd)); % ATM this only makes rasters that show sig control mech response
clIDind = clInfo.id(TblInd);
lngth = length(clIDind);
for a = 1:lngth
    cl = clIDind(a);
    clSel = find(ismember(pclID, cl));
    if chCond == 1
        rasCondSel = find(ismember(consCondNames, MchCondControl) | ismember(consCondNames, MchLasCond));
        label = 'Mech';
    else
        rasCondSel = find(ismember(consCondNames, LasCondControl) | ismember(consCondNames, MchLasCond));
        label = 'Laser';
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
            ax(lidx) = subplot(Nrcond, Nrcl, lidx);
            title(ax(lidx),sprintf('%s cl:%s',rasCondNames{cc},pclID{clSel(ccl)}))
            plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                timeLapse, fs,'',ax(lidx));
            ax(lidx).YAxisLocation = 'origin';ax(lidx).YAxis.TickValues = Nma;
            ax(lidx).YAxis.Label.String = Nma;
            ax(lidx).YAxis.Label.Position =...
                [timeLapse(1)-timeLapse(1)*0.65, Nma,0];
            ax(lidx).XAxis.TickLabels =...
                cellfun(@(x) str2num(x)*1e3, ax(lidx).XAxis.TickLabels,...
                'UniformOutput', 0);
            xlabel(ax(lidx), 'Time [ms]')
            initSub = 0;
            optsRect = {'EdgeColor','none','FaceColor','none'};
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
    rasFigName = ['Unit_', cell2mat(cl), '_', label, '_', num2str(pwr), 'mW_Raster']; 
    configureFigureToPDF (rasFig);
    savefig(rasFig,fullfile(rasterDir, [rasFigName, '.fig']));
end
end
close all
%% ISI Hists
ConsConds = Conditions(3:end);
nCond = length(ConsConds);
% sortedData = sortedData(:,1);
for ChCond = 1:nCond
    ISIhist(ChCond).name = ConsConds(ChCond).name;
    ISIhist(ChCond).Vals(1).name = 'Spontaneous';
    ISIhist(ChCond).Vals(2).name = 'Evoked';
end

for ChCond = 1:nCond
    name = [ConsConds(ChCond).name, '_MR'];
    
    if contains(ConsConds(ChCond).name, 'Laser_Con', 'IgnoreCase', true)
        name = [ConsConds(ChCond).name, '_LR'];
    elseif contains(ConsConds(ChCond).name, ['ch_Laser'], 'IgnoreCase', true)
       name = [ConsConds(ChCond-1).name, '_MR']; % make this more robust to match laser intensity
    end
    
    Ind = find(clInfo.(name));
    ID = clInfo.id(Ind);
    id = find(ismember(gclID, ID))'; % contains will give multiple units when looking for e.g. cl45
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
        [~, isiStack] = getStacks(false,ConsConds(ChCond).Triggers, onOffStr,...
            Window,fs,fs,[],ISIspar(id,:));
        % timelapse becomes spontaneousWindow for pre-trigger, and responseWindow
        % for post
        for histInd = 1: length(id)
            
            figure('Visible','off');
            hisi = histogram(log(isiStack(histInd,:,:)), 'BinEdges', log(1/fs):0.01:log(100));
            ISIhist(ChCond).Vals(wIndex).cts{histInd} = hisi.BinCounts;
            ISIhist(ChCond).Vals(wIndex).bns{histInd} = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
        end
    end
end
%% ISIs and CumISIs
for b = 1:length(ChCond)
    for a = 1:ISIhist(ChCond).Vals(wIndex).cts
        ISIhist(ChCond).Vals(1).ISI{a} = ISIhist(ChCond).Vals(1).cts{a}/sum(ISIhist(ChCond).Vals(1).cts{a});
        ISIhist(ChCond).Vals(2).ISI{a} = ISIhist(ChCond).Vals(2).cts{a}/sum(ISIhist(ChCond).Vals(2).cts{a});
        ISIhist(ChCond).Vals(1).CumISI{a} = ISIhist(ChCond).Vals(1).cts{a}/sum(ISIhist(ChCond).Vals(1).cts{a});
        ISIhist(ChCond).Vals(2).CumISI{a} = ISIhist(ChCond).Vals(2).cts{a}/sum(ISIhist(ChCond).Vals(2).cts{a});
    end
end
%% Saving ISIhist
save(fullfile(dataDir,[expName,'_ISIhist.mat']), 'ISIhist', 'ConsConds', '-v7.3');
%% Saline vs CFA Spontaneous
for model = 1:2
    if model == 1
        Hist = SalISIhist;
    else
        Hist = CfaISIhist;
    end
    
    IsiR = length(Hist(1).Vals(1).cumISI);
    for a = 2: length(Hist)
        IsiR = [IsiR, length(Hist(a).Vals(1).cumISI)];
    end
    r = sum(IsiR);
    Isi = zeros(r, length(Hist(1).Vals(1).cumISI{1}));
    c = 1;
    for a = 1:length(Hist)
        for b = 1:length(Hist(a).Vals(1).cumISI)
            Isi(c,:) = Hist(a).Vals(1).cumISI{b};
            c = c +1;
        end
    end
    % Getting rid of NaNs
    Isi(length(Isi(:,1)) + 1,:) = NaN;
    IndNaN = length(Isi(:,1));
    for nInd = 1:length(Isi(:,1))
        if isnan(Isi(nInd,:)) == true
            IndNaN = [IndNaN; nInd];
        end
    end
    IndNaN = sort(IndNaN, 'descend');
    Isi(IndNaN,:) =[];
    if model == 1
        SpontSalISI = Isi;
        cumSalspont = sum(Isi)/length(Isi(:,1));
    else
        SpontCfaISI = Isi;
        cumCFAspont = sum(Isi)/length(Isi(:,1));
    end
end
figure('Color',[1,1,1], 'Name', 'Saline vs CFA Spontaneous Cumulative Fraction');
% Only if the bin widths are constant!!!
plot(SalISIhist(1).Vals(1).bns{1}, cumSalspont);
hold on
plot(CfaISIhist(1).Vals(1).bns{1}, cumCFAspont);
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([log(0.001), 4.5]);
ylim([0, 1]);
legend('Saline', 'CFA');
title('Spontaneous Activity Cumulative Fractions');
fig = gcf;
ax = gca;
ax.FontSize = 20;
%%
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
%% Stats Test and p-Values on Figure
[txt, star] = findKSsignificance(cumSalspont, cumCFAspont);

annotation(figure(gcf),'textbox',...
    [0.21428125 0.438809261300992 0.0568125 0.029768467475193],'String',star,...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure(gcf),'textbox',...
    [0.673958333333333 0.112328008579419 0.203645833333333 0.232874455900804],...
    'String',{txt},...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
%% Saline & CFA Spont vs Evoked (Mech Control)
for model = 1:2
    if model == 1
        Hist = SalISIhist;
    else
        Hist = CfaISIhist;
    end
    evInd =[];
    for e = 1:length(Hist)
        if contains(Hist(e).name, 'Mech_Control', 'IgnoreCase', true) == true
            evInd = [evInd; e];
        end
    end
    IsiR = length(Hist(evInd(1)).Vals(1).cumISI);
    for a = 2:length(evInd)
        IsiR = [IsiR, length(Hist(evInd(a)).Vals(2).cumISI)];
    end
    r = sum(IsiR);
    Isi = zeros(r, length(Hist(evInd(1)).Vals(2).cumISI{1}));
    c = 1;
    for a = 1: length(evInd)
        for b = 1:length(Hist(evInd(a)).Vals(2).cumISI)
            Isi(c,:) = Hist(evInd(a)).Vals(2).cumISI{b};
            c = c + 1;
        end
    end
    % Getting rid of NaNs
    Isi(length(Isi(:,1)) + 1,:) = NaN;
    IndNaN = length(Isi(:,1));
    for nInd = 1:length(Isi(:,1))
        if isnan(Isi(nInd,:)) == true
            IndNaN = [IndNaN; nInd];
        end
    end
    IndNaN = sort(IndNaN, 'descend');
    Isi(IndNaN,:) =[];
    if model == 1
        EvokedSalISI = Isi;
        cumSalevoked = sum(Isi)/length(Isi(:,1));
    else
        EvokedCfaISI = Isi;
        cumCFAevoked = sum(Isi)/length(Isi(:,1));
    end
end

figure('Color',[1,1,1], 'Name', 'Saline & CFA Spontaneous vs Evoked Cumulative Fraction');
subplot(2,1,1)
for model = 1:2
    subplot(2,1,model)
    if model == 1
        Hist = SalISIhist;
        spont = cumSalspont;
        evoked = cumSalevoked;
    else
        Hist = CfaISIhist;
        spont = cumCFAspont;
        evoked = cumCFAevoked;
    end
    % Only if the bin widths are constant!!!
    plot(Hist(1).Vals(1).bns{1}, spont);
    hold on
    plot(Hist(1).Vals(2).bns{1}, evoked);
    ylabel('Cumulative Fraction');
    xlabel('ISI (msecs)');
    %     xlim([log(0.001), 4.5]);
    %     ylim([0, 1]);
    legend('Spontaneous', 'Evoked');
    fig = gcf;
    ax = gca;
    ax.FontSize = 20;
    if model == 1
        title('Saline');
    else
        title('CFA');
    end
    %%
    ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
    %% More Stats
    [txt, star] = findKSsignificance(spont, evoked);
    annotation(figure(gcf),'textbox',...
        [0.21428125, 0.438809261300992*model 0.0568125 0.029768467475193],'String',star,...
        'FontSize',20,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    annotation(figure(gcf),'textbox',...
        [0.673958333333333, 0.112328008579419*model 0.203645833333333 0.232874455900804],...
        'String',{txt},...
        'FontSize',15,...
        'FitBoxToText','off',...
        'EdgeColor','none');
end

%% CumISI Comparing Conditions
r = 1;
for model = 1:2
    if model == 1
        Histy = SalISIhist;
        modName = 'Saline';
    else
        Histy = CfaISIhist;
        modName = 'CFA';
    end
    for wIndex = 1:2
        if wIndex == 1
            wName = 'Spontaneous';
        else
            wName = 'Evoked';
        end
        
        figure('Color',[1,1,1], 'Name', [modName, ' ', wName, ' Cumulative Fraction']);
         c = 0;
         
        for a = 1:3
            pwr = a*5;
            if a == 1
                chCond = [7 8 9];
            elseif a == 2
                chCond = [1 2 3];
            elseif a == 3
                chCond = [4 5 6];
            end
            subplot(3,1,a)
            hold on
            for chFig = chCond
                Isi = zeros(length(Histy(chFig).Vals(wIndex).cumISI), length(Histy(chFig).Vals(wIndex).cumISI{1}));
                for cInd = 1:length(Histy(chFig).Vals(wIndex).cumISI)
                    Isi(cInd,:) =  Histy(chFig).Vals(wIndex).cumISI{cInd};
                end
                % Getting rid of NaNs
                Isi(length(Histy(chFig).Vals(wIndex).cumISI) + 1,:) = NaN;
                IndNaN = (length(Histy(chFig).Vals(wIndex).cumISI) + 1);
                for nInd = 1:(length(Histy(chFig).Vals(wIndex).cumISI) + 1)
                    if isnan(Isi(nInd,:)) == true
                        IndNaN = [IndNaN; nInd];
                    end
                end
                IndNaN = sort(IndNaN, 'descend');
                Isi(IndNaN,:) =[];
                cumISI = sum(Isi)/length(Histy(chFig).Vals(wIndex).cumISI);
                % Only if the bin widths are constant!!!
                plot(CfaISIhist(1).Vals(1).bns{1}, cumISI)
                if model == 1
                    SalISIhist(chFig).Vals(wIndex).cumsum = cumISI;
                    Histy = SalISIhist;
                elseif model == 2
                    CfaISIhist(chFig).Vals(wIndex).cumsum = cumISI;
                    Histy = CfaISIhist;
                end
            end
            % Results(c).name = kstest2(cumISI(chCond(1)).Condition(wIndex).Vals;
            ylabel('Cumulative Fraction');
            xlabel('ISI (msecs)');
            xlim([log(0.001), 4.5]);
            ylim([0, 1]);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            fig = gcf;
            ax = gca;
            ax.FontSize = 20;
            ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
           
            
            for a = 1:length(chCond)-1
                for b  = a+1: length(chCond)
                    Results(r).name = [modName, ' ', Histy(a + c).name, ' vs ', Histy(b + c).name, ' ', wName];
                    Results(r).kstest2 = findKSsignificance(Histy(a + c).Vals(wIndex).cumsum,Histy(b + c).Vals(wIndex).cumsum);
                    r = r + 1;
                end
            end
            
            
          c = c + 3;  
        end
        
    end
end
for a = 1:length(Results)
fprintf(1, [Results(a).name,' Kolgorov-Smirnov Test ', num2str(Results(a).kstest2), ' \n \n'])
end