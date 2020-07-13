% 30.08.19 Jittering analysis stolen from Emilio by Ross
% clearvars
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

%% User controlling variables RAM Change defaults
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

%% Evening out Condition trials  RAM: comment!
% turn into function
for a = 1: length(Conditions)
    sZ(1,a) = length(Conditions(a).Triggers);
end
[r, c] = min(sZ);
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
Conditions(length(Conditions)).Triggers = concatCond;
clear concatCond;



%this can be put into a function or replaced?
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
save(fullfile(dataDir,[expName,'_Variables.mat']),'consCondNames','Counts', 'Results', 'responseWindow','Triggers', '-v7.3');

%RAM: shove all of this previous stuff into functions and then replace
%with one or more function calls.

%% Getting cluster info and adding variables to table
clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
chanMap = readNPY(fullfile(dataDir,'channel_map.npy'));
chanPos = readNPY(fullfile(dataDir,'channel_positions.npy'));
[m,b] = lineariz(chanPos(:,1), 6, 1);
shank = m*chanPos(:,1) + b;
shank = round(shank);
shMap = containers.Map(chanMap, shank);
setShank = @(x) shMap(x);
clInfo.shank = arrayfun(setShank, clInfo.channel);
tb = size(clInfo);
sz = tb(1);
ActiveUnit = false(sz,1);
clInfo = addvars(clInfo,ActiveUnit,'NewVariableNames','ActiveUnit','After','id');
clInfo{gclID, 'ActiveUnit'} = true;

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
% writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));

%% Determining nShanks

tShanks = sum(clInfo.shank == [1:100]);
nShanks = sum(tShanks ~= false);

%% Plotting spontaneous activity rates
% RAM separate calculations and plotting

% rW = responseWindow(2)-responseWindow(1);
c = 1;
d = 0;
for shankNo = 1:nShanks
    
    index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
    
    for a = 1: length(consCondNames)
        Spont = [consCondNames{1,a}, '_Counts_Spont'];
        SpontaneousBox{d+a} = clInfo.(Spont)(index)/rW;
        SpontMed(1,(d+a)) = median(SpontaneousBox{a});
        SLabels = {consCondNames{1,1}, consCondNames{1,a}};
    end
    
    for a = 2:length(consCondNames)
        figure('Name', ['Spontaneous_Rates_Shank_', num2str(shankNo)], 'Color', 'white');
        boxplot([SpontaneousBox{1,d+1}, SpontaneousBox{1,d+a}]);
        hold on
        for i = 1:length(SpontaneousBox{1,d+1})
            plot([1, 2],[SpontaneousBox{1,d+1}(i,1) SpontaneousBox{1,d+a}(i,1)],'-o', 'color', [0.9,0.9,0.9]);
        end
        SLabels = {consCondNames{1,1}, consCondNames{1,a}};
        title(['Shank ', num2str(shankNo)]);
        ylabel('Firing Rate (Hz)');
        xticklabels(SLabels)
        
        ax = gca;
        ax.FontSize = 14;
        ax = gca;
        ax.FontSize = 14;
    end
    
    % configureFigureToPDF(SpontaneousBox);
    
    
    
    % Getting Wilcoxon rank sums for box plots
    
    for a = 1: length(consCondNames) - 1
        for b = (a + 1): length(consCondNames)
            SpontRS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
            SpontRS(c).SignRank = signrank(SpontaneousBox{d+a}, SpontaneousBox{d+b});
            if SpontRS(c).SignRank <= 0.05
                SpontRS(c).Signifcant = true;
            end
            
            c = c + 1;
        end
    end
    d = d + length(consCondNames);
end
clear SLabels
%% Plotting normalised spontaneous activity rates
% rW = responseWindow(2)-responseWindow(1);
d = 0;

for shankNo = 1:nShanks
    figure('Name', ['Normalised_Spontaneous_Rates_Shank_', num2str(shankNo)], 'Color', 'white');
    for a = 2: length(consCondNames)
        hold on
        plot([0, length(consCondNames)], [0, 0], 'color', [0,0,0]);
        for i = 1:length(SpontaneousBox{1,d+1})
            delta{d+a}(i,1) = ((SpontaneousBox{1,d+1}(i,1)/(SpontaneousBox{1,d+a}(i,1)))-1)*100;
            if delta{d+a}(i,1) > 1
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [0,0,1]);
            elseif delta{d+a}(i,1) < 1
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [1,0,0]);
            else
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [0,0,0]);
            end
            
        end
        
        title(['Normalised Spontaneous Rates: Shank ', num2str(shankNo)]);
        ylabel('% FR Change');
        ax = gca;
        ax.FontSize = 20;
        SLabels{a} = consCondNames{1,a};
        xticklabels(SLabels);
    end
    d = d + length(consCondNames);
end
s
%% making a histogram of delta values for 2?? conditions (check which ones these are)
figure
nBins=200;
[h bins]=hist(delta,nBins);
plot(bins,h,'linewidth', 2)
% configureFigureToPDF(SpontaneousBox);
legend(SLabels(2:end))
legend boxoff

%% small aside for noise correlation population  --should be put into a different script, better labels
close all

Rs=[];Re=[];
for i=1:numel(Counts(:,1))
    C=Counts{i,1};
    %correlation matrix for population
    Rs(:,:,i)=corrcoef(C');
end
for i=1:numel(Counts(:,2))
    C=Counts{i,2};
    %correlation matrix for population
    Re(:,:,i)=corrcoef(C');
end
figure
pcolor([Rs(:,:);Re(:,:)]),shading flat

%[colormap] = cbrewer(ctype, cname, ncol [, interp_method])  FIX THIS
%REBECCA
%plotting stuff

figure
[h bins]=hist([Rs(:) Re(:)],100)
plot(bins,h,'linewidth', 2)
legend('spont', 'evoked')
%
%scatter plot for a neuron pair

n1=1; %neuron 1
n2=2; %neuron 2
figure
plot(C(n1,:),C(n2,:),'.','markersize',15)
xlabel 'n1'
ylabel 'n2'
axis square
grid on




%% Plotting evoked activity rates

% rW = responseWindow(2)-responseWindow(1);
c = 1;
d = 0;
for shankNo = 1:nShanks
    
    index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
    
    for a = 1: length(consCondNames)
        Spont = [consCondNames{1,a}, '_Counts_Evoked'];
        EvokedBox{d+a} = clInfo.(Spont)(index)/rW;
        Evoked_Med(1,(d+a)) = median(EvokedBox{d+a});
        ELabels = {consCondNames{1,1}, consCondNames{1,a}};
    end
    
    for a = 2:length(consCondNames)
        figure('Name', ['Evoked_Rates_Shank_', num2str(shankNo)], 'Color', 'white');
        boxplot([EvokedBox{1,d+1}, EvokedBox{1,d+a}]);
        hold on
        for i = 1:length(EvokedBox{1,d+1})
            plot([1, 2],[EvokedBox{1,d+1}(i,1) EvokedBox{1,d+a}(i,1)],'-o', 'color', [0.9,0.9,0.9]);
        end
        ELabels = {consCondNames{1,1}, consCondNames{1,a}};
        title(['Shank ', num2str(shankNo)]);
        ylabel('Firing Rate (Hz)');
        xticklabels(ELabels)
        
        ax = gca;
        ax.FontSize = 14;
        ax = gca;
        ax.FontSize = 14;
    end
    
    % configureFigureToPDF(SpontaneousBox);
    
    
    
    % Getting Wilcoxon rank sums for box plots
    
    for a = 1: length(consCondNames) - 1
        for b = (a + 1): length(consCondNames)
            EvokedRS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
            EvokedRS(c).SignRank = signrank(EvokedBox{d+a}, EvokedBox{d+b});
            if EvokedRS(c).SignRank <= 0.05
                EvokedRS(c).Signifcant = true;
            end
            
            c = c + 1;
        end
    end
    d = d + length(consCondNames);
end
clear ELabels

%% Plotting normalised evoked activity rates
% rW = responseWindow(2)-responseWindow(1);
d = 0;

for shankNo = 1:nShanks
    figure('Name', ['Normalised_Evoked_Rates_Shank_', num2str(shankNo)], 'Color', 'white');
    for a = 2: length(consCondNames)
        hold on
        plot([0, length(consCondNames)], [0, 0], 'color', [0,0,0]);
        for i = 1:length(EvokedBox{1,d+1})
            delta{d+a}(i,1) = ((EvokedBox{1,d+1}(i,1)/(EvokedBox{1,d+a}(i,1)))-1)*100;
            if delta{d+a}(i,1) > 1
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [0,0,1]);
            elseif delta{d+a}(i,1) < 1
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [1,0,0]);
            else
                plot((a-1),[delta{d+a}(i,1)-1],'-o', 'color', [0,0,0]);
            end
            
        end
        
        title(['Normalised Evoked Rates: Shank ', num2str(shankNo)]);
        ylabel('% FR Change');
        ax = gca;
        ax.FontSize = 20;
        ELabels{a} = consCondNames{1,a};
        xticklabels(ELabels);
    end
    d = d + length(consCondNames);
end
%% Population MRs for each condition.

% Filter for only mechanically responsive clusters?
ansFilt = questdlg('Would you like to filter for significance?','Filter',...
    'Yes','No','Yes');
if strcmp(ansFilt,'Yes')
    % rW = responseWindow(2)-responseWindow(1);
    c = 0;
    d = 0;
    for shankNo = 1:nShanks
        figure('Name', ['Filtered_Mechanical_Responses_Shank_', num2str(shankNo)], 'Color', 'white');
        
        
        for a = 1: length(consCondNames)
            Sig = [consCondNames{1,a}, '_MR'];
            index = find(clInfo.(Sig) & clInfo.shank == shankNo);
            Spont = [consCondNames{1,a}, '_Counts_Spont'];
            Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
            SpBox{d+a} = clInfo.(Spont)(index)/rW;
            EvBox{d+a} = clInfo.(Evoked)(index)/rW;
            SpMed(1,(d+a)) = median(SpBox{d+a});
            EvMed(1,(d+a)) = median(EvBox{d+a});
            MLabels{a,1} = consCondNames{1,a};
            subplot(2,length(consCondNames), a);
            boxplot([SpBox{d+a}, EvBox{d+a}]);
            hold on
            for i = 1:numel(EvBox{d+a})
                plot([1, 2],[SpBox{1,d+a}(i,1) EvBox{1,d+a}(i,1)],'-o', 'color', [0.9,0.9,0.9]) ;
            end
            title(consCondNames{a});
            if a == 1
                ylabel('Firing Rate (Hz)');
            end
            xticklabels({'Spont', 'Evoked'});
            ylim([0 25]);
            ax = gca;
            ax.FontSize = 20;
            subplot(2,length(consCondNames), a + length(consCondNames));
            pie([(sum(clInfo.ActiveUnit & clInfo.shank == shankNo) - length(index)), length(index)]);
            labels = {'Unesponsive','Responsive'};
            legend(labels,'Location','southoutside','Orientation','vertical')
            ax = gca;
            ax.FontSize = 15;
            
            if sum(SpBox{d+a}) ~= false && sum(EvBox{d+a} ~= false)
                c = c + 1;
                MechRS(c).name = [consCondNames{1,a}, '_Spont_vs_Evoked_Shank_', num2str(shankNo)];
                MechRS(c).SignRank = signrank(SpBox{d+a}, EvBox{d+a});
                if MechRS(c).SignRank <= 0.05
                    MechRS(c).Signifcant = true;
                end
                
                
            end
            
        end
        d = d + length(consCondNames);
    end
    
else
    % rW = responseWindow(2) - responseWindow(1);
    c = 1;
    d = 0;
    for a = 1:length(consCondNames)
        figure('Name',['Unfiltered MechResponse_', consCondNames{1,a}], 'Color', 'white')
        
        
        for shankNo = 1:nShanks
            index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
            Spont = [consCondNames{1,a}, '_Counts_Spont'];
            SpBox{d+a} = (clInfo.(Spont)(index))/rW;
            Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
            EvBox{d+a} = (clInfo.(Evoked)(index))/rW;
            subplot(1,nShanks, shankNo);
            boxplot([SpBox{d+a}, EvBox{d+a}]);
            hold on
            for i = 1:numel(EvBox{d+a})
                plot([1 2],[SpBox{d+a}(i,1) EvBox{d+a}(i,1)],'-o', 'color', [0.9,0.9,0.9]) ;
            end
            
            if shankNo == 1
                title(['Shank ', num2str(shankNo)]);
                ylabel('Firing Rate (Hz)');
                xticklabels({'Spont', 'Evoked'});
            else
                title([num2str(shankNo)]);
                xticklabels({[], []});
            end
            ylim([0 20]);
            ax = gca;
            ax.FontSize = 13;
            MechRS(c).name = [consCondNames{1,a}, '_Spont_vs_Evoked_Shank_', num2str(shankNo)];
            MechRS(c).SignRank = signrank(SpBox{d+a}, EvBox{d+a});
            if MechRS(c).SignRank <= 0.05
                MechRS(c).Signifcant = true;
            end
            c = c + 1;
        end
    end
end

%% Relative Responses

% Filter for only mechanically responsive clusters?
ansFilt = questdlg('Would you like to filter for significance?','Filter',...
    'Yes','No','Yes');
if strcmp(ansFilt,'Yes')
   % rW = responseWindow(2)-responseWindow(1);
    c = 1;
    d = 0;
    for shankNo = 1:nShanks
        figure('Name', ['Filtered_Relative_Responses_Shank_', num2str(shankNo)], 'Color', 'white');
        
        
        for a = 1: length(consCondNames)
            Sig = [consCondNames{1,a}, '_MR'];
            index = find(clInfo.(Sig) & clInfo.shank == shankNo);
            Spont = [consCondNames{1,a}, '_Counts_Spont'];
            Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
            RRBox{d+a} = (clInfo.(Evoked)(index))- (clInfo.(Spont)(index))/rW;
            RRMed(1,(d+a)) = median(RRBox{d+a});
            RLabels{a,1} = consCondNames{1,a};
            subplot(1,length(consCondNames), a);
            boxplot(RRBox{d+a});
           if a == 1
                title(['Shank ', num2str(shankNo)]);
                ylabel('Relative Responses (Hz)');
            end
            xticklabels(RLabels{a,1})
            ylim([-10 10]);
            ax = gca;
            ax.FontSize = 20;
        end
        
        
        
        % Getting Wilcoxon rank sums for box plots
        
        for a = 1: length(consCondNames) - 1
            for b = (a + 1): length(consCondNames)
                if sum(RRBox{d+a}) ~= false && sum(RRBox{d+b} ~= false)
                    RR_RS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
                    RR_RS(c).SignRank = signrank(RRBox{d+a}, RRBox{d+b});
                    if RR_RS(c).SignRank <= 0.05
                        RR_RS(c).Signifcant = true;
                    end
                end
                
                c = c + 1;
            end
        end
        d = d + length(consCondNames);
        
    end
else
    
    
    % rW = responseWindow(2) - responseWindow(1);
    c = 0;
    d = 0;
    for shankNo = 1:nShanks
        index = find(clInfo.ActiveUnit & clInfo.shank == shankNo);
        
        for a = 1:length(consCondNames)
            Spont = [consCondNames{1,a}, '_Counts_Spont'];
            Evoked = [consCondNames{1,a}, '_Counts_Evoked'];
            RRBox{d+a} = (clInfo.(Evoked)(index))- (clInfo.(Spont)(index))/rW;
        end
        for a = 2: length(consCondNames)
            figure('Name',['Unfiltered_Relative_Responses_Shank_', num2str(shankNo),' ', num2str(consCondNames{1,a})], 'Color', 'white')
            boxplot([RRBox{1,d+1}, RRBox{1,d+a}])
            hold on
            for i = 1:numel(RRBox{1,d+1})
                plot([1 2],[RRBox{1,d+1}(i,1) RRBox{1,d+a}(i,1)],'-o', 'color', [0.9,0.9,0.9]) ;
            end
            title(['Shank ', num2str(shankNo)]);
            ylabel('Firing Rate (Hz)');
            RLabels = {consCondNames{1,1}, consCondNames{1,a}};
            xticklabels(RLabels);
            ylim([-5 20]);
            ax = gca;
            ax.FontSize = 13;
        end
        for a = 1: length(consCondNames) - 1
            for b = (a + 1): length(consCondNames)
                if sum(RRBox{d+a}) ~= false && sum(RRBox{d+b} ~= false)
                    RR_RS(c).name = [consCondNames{1,a},'_vs_', consCondNames{1,b}, '_Shank_', num2str(shankNo)];
                    RR_RS(c).SignRank = signrank(RRBox{d+a}, RRBox{d+b});
                    if RR_RS(c).SignRank <= 0.05
                        RR_RS(c).Signifcant = true;
                    end
                end
                
                c = c + 1;
            end
        end
        d = d + length(consCondNames);
    end
end

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
%% Plotting the population activity

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
goodsIdx = ~badsIdx';
csNames = fieldnames(Triggers);
for ccond = 1:Nccond
    figFileName = sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, Conditions(consideredConditions(ccond)).name, timeLapse*1e3,...
        binSz*1e3, responseWindow*1e3, spontaneousWindow*1e3, onOffStr,...
        orderedStr, filtStr);
    [PSTH, trig, sweeps] = getPSTH(discStack(filterIdx,:,:),timeLapse,...
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
    figs = plotClusterReactivity(PSTH(ordSubs,:),trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};...
        pclID(ordSubs)],...
        strrep(expName,'_','\_'),...
        stims, csNames);
    configureFigureToPDF(figs);
    figs.Children(end).YLabel.String = [figs.Children(end).YLabel.String,...
        sprintf('^{%s}',orderedStr)];
    if ~exist([figFileName,'.pdf'], 'file') || ~exist([figFileName,'.emf'], 'file')
        print(figs,fullfile(figureDir,[figFileName, '.pdf']),'-dpdf','-fillpage')
        print(figs,fullfile(figureDir,[figFileName, '.emf']),'-dmeta')
    end
end



%% rasters for testing
figure
lastlevel=0;
thisCondition=3;
for n=40:50 %numel(clusters)
    if mod(n,2)==0, col='b';else col='k';end
    [R lastlevel]=manyRasters(relativeSpikeTimes(n,thisCondition).sp,col,.9,ppms,0+lastlevel);
end
title(relativeSpikeTimes(1,thisCondition).condition)
grid on

%% to-do
%all files to run Emilio code on SDS for sample file Z:\Ross\Experiments\10mW_CFA_VPL
% compare rasters created by different triggering methods
% compare rasters for problematic PSTHs
% trial-by-trial rasters (grouped by trial rather than neuron)
% can we find actual mechanical stimulus?
% add triggered stimulus (especially light) to plots to check for artifacts
% or units directly driven by light

