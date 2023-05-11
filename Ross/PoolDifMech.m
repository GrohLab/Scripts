%% PoolDifMech

clear
close all
clc
%%
dataDirs = {
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\1.8.22\KS3__Nblocks2__9_9__0pt9__20\'
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\4.8.22\KS3__Nblocks2__9_9__0pt9__20\'
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\8.8.22\DifMech\VPL_E1_KS3__Nblocks2__9_9__0pt9__20\'
    };

selUnitDirs = {
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\1.8.22\KS3__Nblocks2__9_9__0pt9__20\m51_ECE_Processing_-10-to-10\PopulationAnalysis\SelectedUnitData.mat'
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\4.8.22\KS3__Nblocks2__9_9__0pt9__20\m50_ECE_Processing_-10-to-10\PopulationAnalysis\SelectedUnitData.mat'
    '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\8.8.22\DifMech\VPL_E1_KS3__Nblocks2__9_9__0pt9__20\m53_VPL_E1_DifMech_ECE_Processing_-10-to-10\PopulationAnalysis\SelectedUnitData.mat'
    };


%% variables to be collected

discStacks = [];
csts = [];
Data = [];

%% collecting trial numbers early
chCond = 1;
nExpts = size(dataDirs, 1);
trialnumbers = [];
for cexpt = 1:nExpts
    dataDir = dataDirs{cexpt};
    Conditions = dir([dataDir, '*analysis*.mat']);
    [r, ~] = size(Conditions);
    if r ~= 1
        fprintf(['multiple or no analysis files in directory...\n' ...
            'choose one \n'])
        return
    end
    load([dataDir, Conditions.name]);
    trialnumbers = [trialnumbers, size(Conditions(1).Triggers,chCond)];
end

constrials = min(trialnumbers);
%% User controlling variables
% Time lapse, bin size, and spontaneous and response windows
promptStrings = {'Viewing window (time lapse) [s]:','Response window [s]',...
    'Bin size [s]:'};
defInputs = {'-2, 6', '0.1, 5', '0.01'};
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




%% loading variables for a given experiment - Starting the for-loop

for cexpt = 1:nExpts
    dataDir = dataDirs{cexpt};
    selUnits = selUnitDirs{cexpt};

    Conditions = dir([dataDir, '*analysis*.mat']);
    [r, ~] = size(Conditions);
    if r ~= 1
        fprintf(['multiple or no analysis files in directory...\n' ...
            'choose one \n'])
        return
    end
    load([dataDir, Conditions.name]);

    % homogenising trial numbers
    Conditions(chCond).Triggers = Conditions(chCond).Triggers(1:constrials,:);


    sortedData = dir([dataDir, '*all_channels*.mat']);
    [r, ~] = size(sortedData);
    if r ~= 1
        fprintf(['multiple or no all_channels files in directory...\n' ...
            'choose one \n'])
        return
    end
    load([dataDir, sortedData.name]);

    clInfo = getClusterInfo([dataDir filesep 'cluster_info.tsv']);

    expID =  dir([dataDir, '*expParams*.mat']);
    [r, ~] = size(expID);
    if r ~= 1
        fprintf(['multiple or no expParams files in directory...\n'])
        return
    end
    load([dataDir, expID.name], 'expID');

    load(selUnits, 'unit_ids')


    %% renaming units to avoid confusion


    for i = 1:length(sortedData)
        sortedData{i,1} = [expID, '_unit_', sortedData{i,1}];
    end

    for i = 1:height(clInfo)
        clInfo.id{i,1} = [expID, '_unit_', clInfo.id{i,1}];
        clInfo.Properties.RowNames{i,1} = clInfo.id{i,1};
    end

    for i = 1:length(unit_ids)
        unit_ids{i} = [expID, '_unit_', unit_ids{i,1}];
    end



    %% global variables

    goods = find(ismember(sortedData(:,1), unit_ids));

    Triggers.MechStim = Triggers.MechStim * -1;

    % Number of total samples
    Ns = min(structfun(@numel,Triggers));
    % Total duration of the recording
    Nt = Ns/fs;

    gclID = sortedData(goods,1);
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



    Data = [Data; clInfo(:,1:13)];

    %% Constructing the stack out of the user's choice

   
    onOffStr = 'on';

    % discStack - dicrete stack has a logical nature
    % cst - continuous stack has a numerical nature
    % Both of these stacks have the same number of time samples and trigger
    % points. They differ only in the number of considered events.

    [discStack, cst] = getStacks(spkLog,Conditions(chCond).Triggers,onOffStr,...
        timeLapse,fs,fs,spkSubs,continuousSignals);
    nTrials = size(discStack, 3);
    if cexpt == 1
        nTrialsStacks = size(discStack, 3);
    end


    tic
    while nTrials > nTrialsStacks
        discStack(:,:,end) = [];
    end
    toc


    tic
    while nTrials < nTrialsStacks
        discStacks(:,:,end) = [];
    end
    toc

    discStacks = cat(1, discStacks, discStack);
    csts = cat(3, csts, cst);
end

discStack = discStacks; clear discStacks
[Ne, Nt, NTa] = size(discStacks);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);


consideredConditions = find(~ismember(1:length(Conditions), chCond));

%% Boolean flags
delayFlags = false(NTa,Nccond);
counter2 = 1;
for ccond = consideredConditions
    delayFlags(:,counter2) = ismember(Conditions(chCond).Triggers(:,1),...
        Conditions(ccond).Triggers(:,1));
    counter2 = counter2 + 1;
end
Na = sum(delayFlags,1);