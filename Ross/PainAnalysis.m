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
    gclID = sortedData(goods,1);
    badsIdx = StepWaveform.subs2idx(bads,size(sortedData,1));
    % Logical spike trace for the first good cluster
    spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % Subscript column vectors for the rest good clusters
    spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods(2:end),2),...
        'UniformOutput',false);
    % Number of good clusters 
    Ncl = numel(goods);
    
   % continuousSignals = 
    
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

%% Evening out Condition trials

for a = 1: length(Conditions)
    sZ(1,a) = length(Conditions(a).Triggers);
end
[r, c] = min(sZ);
minVal = min(sZ);
for a = 1: length(Conditions)
    if a == c
    a = a + 1;
    else
        Conditions(a).Difference = length(Conditions(a).Triggers) - length(Conditions(c).Triggers);
        Ind = randperm(length(Conditions(a).Triggers), Conditions(a).Difference);
        Conditions(a).Triggers(Ind',:) = [];
    end
end
Conditions(length(Conditions)).name = 'All';
concatCond = Conditions(1).Triggers;
for a = 2:length(Conditions) - 1
    concatCond = [concatCond; Conditions(a).Triggers];
end

Conditions(length(Conditions)).Triggers = concatCond; clear concatCond;




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
[~, isiStack] = getStacks(spkLog,Conditions(chCond).Triggers, onOffStr,...
    timeLapse,fs,fs,[],ISIspar);
% [dst, cst] = getStacks(spkLog, allWhiskersPlusLaserControl,...
%     'on',timeLapse,fs,fs,[spkSubs;{Conditions(allLaserStimulus).Triggers}],...
%     continuousSignals);
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
figure; imagesc(delayFlags);
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
figs = scatterSignificance(Results, Counts, consCondNames, delta_t, sortedData(goods,1));
configureFigureToPDF(figs);
stFigBasename = fullfile(figureDir,[expName,' ']);
stFigSubfix = sprintf(' Stat RW%.1f-%.1fms',responseWindow(1)*1e3,...
    responseWindow(2)*1e3);
ccn = 1;
for cc = indCondSubs
    stFigName = [stFigBasename, consCondNames{ccn}, stFigSubfix];
    ccn = ccn + 1;
    if ~exist([stFigName,'.*'],'file')
        print(figs(cc),[stFigName,'.pdf'],'-dpdf','-fillpage')
        print(figs(cc),[stFigName,'.emf'],'-dmeta')
    end
end

%% Saving variables
save(fullfile(dataDir,[expName,'_Variables.mat']),'consCondNames','Counts', 'Results', 'responseWindow','-v7.3');

%% Getting cluster info and adding variables to table
clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
% clInfo = addvars(clInfo,~badsIdx','NewVariableNames','ActiveUnit','After','id');
% clInfo.shank = arrayfun(setShank, clInfo.channel);


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

sZ = size(Counts);
     d = length(consCondNames);
     e = 1;
     f = length(consCondNames);
for a = 1:length(consCondNames) - 1
    for b = (a + 1): length(consCondNames)
        
        for c = 1: sZ(1,2)             
        clInfo{clInfo.ActiveUnit == true,[consCondNames{1,a},'_vs_', consCondNames{1,b}, '_', Results(e).Activity(c).Type, '_Response']} = Results(e).Activity(c).Pvalues < 0.05;
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

%% Accessing different units from table

% IDs for active units
for a = 1: length(consCondNames)
    IdActive.Clusters = find(clInfo.ActiveUnit);
end

% IDs for mechanically responsive clusters per condtion
for a = 1: length(consCondNames)
    IdMR(a).name = [consCondNames{1,a}, '_MR']; IdMR(a).Clusters = find(clInfo.(IdMR(a).name));
end

% IDs for between condition sig differences in spontaneous and evoked FRs
d = 0;
for a = 1:length(consCondNames) - 1
    for b = (a + 1): length(consCondNames)
        
       for c = 1:2
           IdVs(d+c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_', Results(1).Activity(c).Type, '_Response'];
           IdVs(d+c).Clusters = find(clInfo.(IdVs(d+c).name));
       end
      d = d + length(consCondNames);
    end
end

%% Determining nShanks

tShanks = sum(clInfo.shank == [1:100]);
nShanks = sum(tShanks ~= false);

%% Plotting spontaneous activity rates

rW = responseWindow(2)-responseWindow(1);
c = 1;
d = 0;
figure('Name', 'Spontaneous_Rates', 'Color', 'white');
for shankNo = 1:nShanks
    
index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
    
    for a = 1: length(consCondNames)

        Spont = [consCondNames{1,a}, '_Counts_Spont'];
        SpontaneousBox(:,a) = clInfo.(Spont)(index);
        SpontMed(1,(d+a)) = median(SpontaneousBox(:,a))/rW;
        SLabels{a,1} = consCondNames{1,a};
    end
    SpontaneousBox = SpontaneousBox/rW; 
    
    subplot(1, nShanks, shankNo);
    boxplot(SpontaneousBox);
    if shankNo == 1
                title(['Shank ', num2str(shankNo)]);
                ylabel('Firing Rate (Hz)');
                xticklabels(SLabels)
            else
                title([num2str(shankNo)]);
                xticklabels({[], []});
    end
    ylim([0 25]);
    ax = gca; 
    ax.FontSize = 12;
    ax = gca; 
    ax.FontSize = 12;
    % configureFigureToPDF(SpontaneousBox);


    % Getting Wilcoxon rank sums for box plots
    
    for a = 1: length(consCondNames) - 1
        for b = (a + 1): length(consCondNames)        
            SpontRS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
            SpontRS(c).RankSum = ranksum(SpontaneousBox(:,a), SpontaneousBox(:,b));
            if SpontRS(c).RankSum <= 0.05
                SpontRS(c).Signifcant = true;
            end

            c = c + 1;
        end
    end 
    d = d + length(consCondNames);
    clear SpontaneousBox
end

%% Plotting evoked activity rates

rW = responseWindow(2)-responseWindow(1);
c = 1;
d = 0;
figure('Name', 'Evoked_Rates', 'Color', 'white');
for shankNo = 1:nShanks
    
index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
    
    for a = 1: length(consCondNames)

        Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
        EvokedBox(:,a) = clInfo.(Evoked)(index);
        EvokedMed(1,(d+a)) = median(EvokedBox(:,a))/rW;
        ELabels{a,1} = consCondNames{1,a};
    end
    EvokedBox = EvokedBox/rW; 
    
    subplot(1, nShanks, shankNo);
    boxplot(EvokedBox);
    if shankNo == 1
                title(['Shank ', num2str(shankNo)]);
                ylabel('Firing Rate (Hz)');
                xticklabels(ELabels)
            else
                title([num2str(shankNo)]);
                xticklabels({[], []});
    end
    ylim([0 25]);
    ax = gca; 
    ax.FontSize = 12;
    % configureFigureToPDF(EvokedBox);


    % Getting Wilcoxon rank sums for box plots
    
    for a = 1: length(consCondNames) - 1
        for b = (a + 1): length(consCondNames)        
            EvRS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
            EvRS(c).RankSum = ranksum(EvokedBox(:,a), EvokedBox(:,b));
            if EvRS(c).RankSum <= 0.05
                EvRS(c).Signifcant = true;
            end

            c = c + 1;
        end
    end 
    d = d + length(consCondNames);
    clear EvokedBox
end


%% Population MRs for each condition (unfiltered for mech significance).

rW = responseWindow(2)-responseWindow(1);
c = 0;
for a = 1:length(consCondNames)
    figure('Name',['MechResponse_', consCondNames{1,a}], 'Color', 'white') 
      

        for shankNo = 1:nShanks
            index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
            c = c + 1;
            Spont = [consCondNames{1,a}, '_Counts_Spont'];
            SpBox = clInfo.(Spont)(index);
            Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
            EvBox = clInfo.(Evoked)(index);
            subplot(1,nShanks, shankNo);
            boxplot([SpBox, EvBox]);
            
            if shankNo == 1
                title(['Shank ', num2str(shankNo)]);
                ylabel('Firing Rate (Hz)');
                xticklabels({'Spont', 'Evoked'});
            else
                title([num2str(shankNo)]);
                xticklabels({[], []});
            end
            ylim([0 25]);
            ax = gca; 
            ax.FontSize = 12;
            MechRS(c).name = [consCondNames{1,a}, '_Spont_vs_Evoked_Shank_', num2str(shankNo)];
            MechRS(c).RankSum = ranksum(SpBox, EvBox);
            if MechRS(c).RankSum <= 0.05
                    MechRS(c).Signifcant = true;
            end

            clear SpBox; clear EvBox;
        end
end

    
    