% PAIN
% 30.08.19 Jittering analysis by using the Data Explorer.
%clearvars
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

Triggers.MechStim = Triggers.MechStim * -1;

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
% clInfo = getClusterInfo([dataDir, '\cluster_info.tsv']);
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
% try
%     [~, isiStack] = getStacks(spkLog,Conditions(chCond).Triggers, onOffStr,...
%         timeLapse,fs,fs,[],ISIspar);
% catch
%     fprintf(1, 'Perhaps there''s not enough memory to cope with this ISIs\n');
%     % General ISI
%
% end
% %[dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
% %    'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
% %    continuousSignals);
% if ~exist(isiFile,'file')
%     fprintf(1,'Saving the inter-spike intervals for each cluster... ');
%     save(isiFile,'ISIspar','isiStack','-v7.3')
%     fprintf(1,'Done!\n')
% end
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
    savefig(Figs(cc),fullfile(figureDir, [altCondNames, stFigSubfix '.fig']));
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

fprintf('%d responding clusters:\n', Nwru);
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
%% Add modulation index

ml = clInfo.Mech_Laser_5mW_Rate_Evoked;
mc = clInfo.Mech_Control_5mW_Rate_Evoked;
modInd = (ml - mc)./(ml + mc);

if ~any(ismember(clInfo.Properties.VariableNames,'modIndex'))
    clInfo = addvars(clInfo,modInd,'After','shank',...
        'NewVariableNames','modIndex');
end
% writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));
%% adding abs_depth to table

promptStrings = {'How deep was your probe? (um)'};
defInputs = {'1400'};
dpth = inputdlg(promptStrings,'Inputs', [1, 30],defInputs);

dpth = str2double(dpth{1});


abs_depth = dpth - clInfo.depth;
% dif = abs(dpth) - 1275;
% abs_depth = clInfo.depth + dif;




clInfo = addvars(clInfo,abs_depth,'After','depth',...
    'NewVariableNames','abs_depth');

writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));

%% Depth vs Modulation Plots
clInd = ismember(clInfo.id, gclID);
depths = table(clInfo.id(clInd), clInfo.abs_depth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
depth = depths.Var2;
jitter = randi(20,size(depth));
depth = depth + jitter; % adding small jitter for visualisation
depth = -1*depth;
for chc = 1:length(consideredConditions)
    spontCounts = Counts{chc,1};
    evCounts = Counts{chc,2};
    c = size(spontCounts);
    c = c(2);
    ordSpontCounts = zeros(length(ID), c);
    ordEvCounts = zeros(length(ID), c);
    for cl = 1:length(ID)
        ind = find(ismember(gclID, ID(cl)));
        ordSpontCounts(cl,:) = spontCounts(ind,:);
        ordEvCounts(cl,:) = evCounts(ind,:);
    end
    spWindow = spontaneousWindow(2)-spontaneousWindow(1);
    evWindow = responseWindow(2)-responseWindow(1);
    medSpontRate = median(ordSpontCounts,2)/spWindow;
    medEvokedRate = median(ordEvCounts,2)/evWindow;
    deltaRate = zeros(length(depth), 2);
    deltaRate(:,2) = medEvokedRate - medSpontRate;
    normDeltaRate = zeros(length(depth), 2);
    f = (medEvokedRate + medSpontRate);
    normDeltaRate(:,2) = (medEvokedRate - medSpontRate)./f;
    increased{chc} = ID(deltaRate(:,2)>0);
    decreased{chc} = ID(deltaRate(:,2)<0);
    fig = figure('Name', 'Depth vs Modulation', 'Color', 'White');
    ax = gca;
    hold on
    for cl = 1:length(ID)
        dr = deltaRate(cl,:);
        dpth = [depth(cl), depth(cl)];
        plot(dr, dpth, 'LineStyle', '-', 'Marker', 'none', 'Color', [0, 0, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5]);
        
    end
    
    
    %     legend on
    %plot(deltaRate(:,2), depth, 'LineStyle', '-', 'Marker', 'none', 'Color', [0.5, 1, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5])
    title(['Depth vs Change in Firing Rate: ', consCondNames{chc}], 'Interpreter', 'none')
    ylabel( 'Depth_{(\mum)}', 'Interpreter', 'tex');
    xlabel('\Delta Firing Rate_{ Hz}', 'Interpreter', 'tex');
    %ax.XLim = [-1.1*max(deltaRate(:,2)), 1.1*max(deltaRate(:,2))];
    ax.XLim = [-50, 50];
    
    hold on
    plot(deltaRate(:,1), linspace(ax.YLim(1),ax.YLim(2), length(deltaRate)), 'LineStyle', '-', 'Color', 'k');
    ax.FontSize = 20;
    %     ax.Legend.String = {'Late Responders'};
    %     savefig(fig, fullfile('Z:\Ross\Experiments\smrxFiles\16.12.20\12.6.21\Figures\', [consCondNames{chc}, '_LateResponders_Depth_vs_Modulation.fig']));
end

%% Get significantly different clusters
gcans = questdlg(['Do you want to get the waveforms of the good clusters?'], 'Waveforms', 'Yes', 'No', 'No');
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
            'Mechanical',reshape(mean(cst((1),:,delayFlags(:,cdel)),3),...
            1,Nt),...
            'MechPressure',reshape(mean(cst((2),:,delayFlags(:,cdel)),3),...
            1,Nt),...
            'Laser',reshape(mean(cst(lrFlag,:,delayFlags(:,cdel)),3),...
            1,Nt),'TimeAxis',(0:Nt-1)/fs + timeLapse(1));
        cdel = cdel + 1;
    end
    %save(fullfile(dataDir,[expName,'_analysis.mat']),'Conditions','-append')
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
    %     if ~existFlag
    %         fID = fopen(csvFileName,'w');
    %         fprintf(fID,'%s, %s\n','Cluster ID','Relative spike times [ms]');
    %     end
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
    %     if ~exist([figFileName,'.pdf'], 'file')
    %         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.pdf']),...
    %             '-dpdf','-fillpage')
    %     end
    %     if ~exist([figFileName,'.emf'], 'file')
    %         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.emf']),...
    %             '-dmeta')
    %     end
    
end
% for a = 1:length(consideredConditions)
%     savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH_0.001binSz.fig']));
% end


%% depthPSTH

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
csNames = fieldnames(Triggers);
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
    %     psthFigs(ccond).Children(end).YLabel.String =...
    %         [psthFigs(ccond).Children(end).YLabel.String,...
    %         sprintf('^{%s}',orderedStr)];
    %     if ~exist([figFileName,'.pdf'], 'file')
    %         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.pdf']),...
    %             '-dpdf','-fillpage')
    %     end
    %     if ~exist([figFileName,'.emf'], 'file')
    %         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.emf']),...
    %             '-dmeta')
    %     end
    
end
% for a = 1:length(consideredConditions)
%     savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH_0.001binSz.fig']));
% end



%% PopPSTH Comparisons - By Condition


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
    %     ax1 = subplot(totlX,1,1:2,'Parent',fig);
    
    clrDivider = nFilters-1;
    if nFilters <= 1
        clrDivider = nFilters;
    end
    %datclr =  [0:0.8/(nFilters-1):0.8]' * ones(1, 3);
    datclr = [green; red; blue];    %[1,0,0; 0,1,0; 0,0,0; 0,0.5,1];
    
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
        %ax(a).YLim = [0, 50];
        %ax(a).YLabel.String = 'Firing Rate [Hz]';
        ax(a).Box = 'off';
        ax(a).ClippingStyle = 'rectangle';
        legend(ax(a),'show','Location','northeast')
        ruNames = [{'L4'}, {'L5'}, {'L6'}];
        ax(a).Legend.String = ruNames{a};
        ax(a).YAxis(1).Label.String = axLabel;
        
    end
    ax(1).Title.String = consCondNames(ccond);
    %ax1.Title.String = consCondNames(ccond);
    % %if fthAxFlag
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
        
        %                    ax2.Children(1).Color = defineColorForStimuli(IDs(cs));
        
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
    
    
    
    
    
    
    
    
    
    
    %      [Ncl, Npt] = size(PSTH);
    %
    %     psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
    %     trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
    %     % clr = defineColorForStimuli(IDe);
    %     popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps);
    %     plot(ax1,psthTX,popPSTH,'DisplayName','Population PSTH')
    %     axLabel = 'Population activity';
    %      yyaxis(ax1,'right')
    %     plot(ax1,trigTX,trig,'LineWidth',1.5,...
    %         'LineStyle',':')%,'Color',clr,'DisplayName',IDe{1},...
    %
    % Formatting the population PSTH plot
    
    ax2.YAxis(1).Limits = [0,1.01];
    ax2.Legend.String = csNames;
    %
    %
    
    %     title(ax1,[expName,sprintf(' %d trials',sweeps)], 'Interpreter', 'none')
    %
    %
    clear PSTH
    clear PSTHarray
    
    
end
% end



%% PopPSTH By Group


red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];
purple = [0.5,0,0.5];
colour = [red; green; blue; purple];


transparency = [1, 0.25]; %1:-1/length(consideredConditions):0.5;

csNames = fieldnames(Triggers);
% ruIdx = filterIdx;
axLabel = 'Firing Rate [Hz]';
ruIdxDim = size(ruIdx);
nFilters = ruIdxDim(2);
for grp = 1: nFilters

    
    fig = figure('Name',[expName,'_', ruNames{grp}],'Color',[1,1,1]);

    hold on
    fthAxFlag = false;
    if ~exist('stims','var')
        totlX = 2;
        
    elseif ~isempty(stims)
        fthAxFlag = true;
        totlX = 3;
    end
    ax1 = subplot(totlX,1,1:2,'Parent',fig);
    
%     
%     datclr =  [0:0.8/(nFilters-1):0.8]' * ones(1, nFilters);
    
    for ccond = length(consideredConditions):-1:1
        [PSTHarray{grp}(:,:,ccond), trig, sweeps] = getPSTH(discStack(ruIdx(:,grp),:,:),timeLapse,...
            ~delayFlags(:,ccond),binSz,fs);
        
        PSTH = PSTHarray{grp}(:,:,ccond);
        [Ncl, Npt] = size(PSTH);
        PSTHn = PSTH./max(PSTH,[],2);
        
        
        psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
        trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
        popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps * binSz);
        plot(psthTX,popPSTH, 'Color', colour(grp,:))
        hold on
        
    end
    leg = legend;
    leg.String = flip(consCondNames);
    leg.Box = 'off';
    leg.FontName = 'Arial';
    leg.FontSize = 10;
    leg.Location = 'north';
    
    ax1.Title.String = ruNames{grp};
    ax1.YLim = [0, 50];
    % %if fthAxFlag
    ax2 = subplot(totlX,1,3,'Parent',fig);
    
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
        
        %                    ax2.Children(1).Color = defineColorForStimuli(IDs(cs));
        
        if cs == 1
            ax2.NextPlot = 'add';
        end
    end
    legend(ax2,'show','Location','northwest')
    
    ax2.Box = 'off';
    ax2.XLim = [timeLapse(1), timeLapse(2)];
    ax2.XAxis.Visible = 'off';
    ax2.YAxis.Visible = 'off';
    linkaxes([ax1,ax2],'x')
    %end
    
    
    
    
    
    
    
    
    
    
    %      [Ncl, Npt] = size(PSTH);
    %
    %     psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
    %     trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
    %     % clr = defineColorForStimuli(IDe);
    %     popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps);
    %     plot(ax1,psthTX,popPSTH,'DisplayName','Population PSTH')
    %     axLabel = 'Population activity';
    %      yyaxis(ax1,'right')
    %     plot(ax1,trigTX,trig,'LineWidth',1.5,...
    %         'LineStyle',':')%,'Color',clr,'DisplayName',IDe{1},...
    %
    % Formatting the population PSTH plot
    ax1.YAxis(1).Label.String = axLabel;
    ax2.YAxis(1).Limits = [0,1.01];
    ax2.Legend.String = csNames;
    ax = gca;
    leg = legend;
    leg.Box = 'off';
    leg.FontName = 'Arial';
    leg.FontSize = 10;
    %
    %
    ax1.XLabel.String = sprintf('Time_{%.2f ms} [s]',binSz*1e3);
    ax1.XLim = [timeLapse(1), timeLapse(2)];
    ax1.Box = 'off';
    ax1.ClippingStyle = 'rectangle';
    %
    legend(ax1,'show','Location','northeast')
%     ax1.Legend.String = ruNames;
    %     title(ax1,[expName,sprintf(' %d trials',sweeps)], 'Interpreter', 'none')
    %
    %
    clear PSTH
    clear PSTHarray
    
    
end
% end

%% Rasters
% DE_Jittering needs to be unfiltered for significance for this to work!
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
    % MchTblInd = ['Mech_Control_', num2str(pwr), 'mW_MR'];
    % LasTblInd =  [20,21,22];%['Laser_Control_', num2str(pwr), 'mW_LR'];
    % MchCondControl = ['Mech_Control_', num2str(pwr), 'mW'];
    % LasCondControl = ['Laser_Control_', num2str(pwr), 'mW'];
    % MchLasCond = ['Mech_Laser_', num2str(pwr), 'mW'];
    % EffectTblInd = ['Mech_Control_', num2str(pwr), 'mW_vs_Mech_Laser_', num2str(pwr), 'mW_Evoked_Response'];
    % TblInd = find(clInfo.ActiveUnit); % ATM this only makes rasters that show sig control mech response
    % clIDind = clInfo.id(TblInd);
    pwrInd = Power == pwr;
    clIDind =  pclID;
    lngth = length(clIDind);
    for a = 1:lngth
        rng('default');
        cl = clIDind(a);
        clSel = find(ismember(pclID, cl));
        %     if chCond == 1
        %         rasCondSel = find(ismember(consCondNames, MchCondControl) | ismember(consCondNames, MchLasCond));
        %         label = 'Mech';
        %     else
        %         rasCondSel = find(ismember(consCondNames, LasCondControl) | ismember(consCondNames, MchLasCond));
        %         label = 'Laser';
        %     end
        rasCondSel = [1 2 3];;
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
%                 ax(lidx).YAxis.Label.Position =...
%                     [timeLapse(1)-timeLapse(1)*0.65, Nma,0];
                %             ax(lidx).XAxis.TickLabels =...
                %                 cellfun(@(x) (x)*1e3, ax(lidx).XAxis.TickValues,...
                %                 'UniformOutput', 0);
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
                
%                 stmClr =[ 1 0 0 0.25; 0 1 1 0.25];
                
                
                
%                 yyaxis('right');
%                 hold on
%                 for cs = 1:min(r,c)
%                     if exist('IDs','var')
%                         plot(trigTX,stims(:,cs),'LineStyle','-','LineWidth', 1,...
%                             'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
%                     else
%                         plot(trigTX,stims(:,cs),'LineStyle','-','LineWidth',1, 'Color', stmClr(cs,:))
%                     end
                    
                    
                    %                    ax2.Children(1).Color = defineColorForStimuli(IDs(cs));
                    
%                     if cs == 1
%                         ax2.NextPlot = 'add';
%                     end

%                 end
                hold off
                ax = gca;
                
               
%                 ax.YAxis(2).Limits = [0.015, 1];
%                 ax.YAxis(2).Visible = 'off';
                ax.FontName ='Arial';
                ax.FontSize = 25;
                
%                 f=get(gca,'Children');
%                 legend(f)
%                 
                
                
                
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
close all



