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

figureDir = '\\lsdf02.urz.uni-heidelberg.de\sd19b001\PainData\Corrected_Channel_Map\VPL\DifMechPooled';
%% variables to be collected

discStacks = false(0);
csts = [];
Data = [];
sortedDatapooled = [];
goodspooled = [];
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
defInputs = {'-2, 6', '0.1, 3', '0.05'};
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

    if any(ismember(clInfo.Properties.VariableNames,'ActiveUnit'))
        clInfo = removevars(clInfo, 'ActiveUnit');
    end


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
    spkSubs = cellfun(@(x) round(x.*fs),sortedData(goods,2),...
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



    Data = [Data; clInfo(:,1:12)];

    %% Constructing the stack out of the user's choice


    onOffStr = 'on';

    % discStack - dicrete stack has a logical nature
    % cst - continuous stack has a numerical nature
    % Both of these stacks have the same number of time samples and trigger
    % points. They differ only in the number of considered events.

    [~, cst] = getStacks(false,Conditions(chCond).Triggers,onOffStr,...
        timeLapse,fs,fs,spkSubs,continuousSignals);
    if cexpt == 1
        nTriggerscst = size(cst, 1);
    end



    while size(cst, 1) > nTriggerscst
        cst(end,:,:) = [];
    end




    while size(cst, 1) < nTriggerscst
        csts(end,:,:) = [];
    end

    discStack = getStacks(false,Conditions(chCond).Triggers,onOffStr,...
        timeLapse,fs,fs,spkSubs,continuousSignals);
    if cexpt == 1
        discStack(2,:,:) = [];
    else
        discStack(1:2,:,:) = [];
    end

    %% concatenating important variables across expts

    tic
    discStacks = cat(1, discStacks, discStack);
    toc
    csts = cat(1, csts, cst);


    goodsoffset = size(sortedDatapooled,1);
    goodspooled = cat(1, goodspooled, goods + goodsoffset);

    sortedDatapooled = cat(1, sortedDatapooled, sortedData);

end
%%
discStack = discStacks; clear discStacks
cst = [mean(csts([1,3,5],:,:)); mean(csts([2,4,6],:,:))]; clear csts
sortedData = sortedDatapooled; clear sortedDatapooled
goods = goodspooled; clear goodspooled
[Ne, Nt, NTa] = size(discStack);
Ncl = numel(goods);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);

condNames = arrayfun(@(x) x.name,Conditions,'UniformOutput',false);
consideredConditions = find(~ismember(1:length(Conditions), chCond));
Nccond = length(consideredConditions);
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
gclID = sortedData(goods,1);
expName = 'pooled';
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
            1,Nt),'TimeAxis',(0:Nt-1)/fs + timeLapse(1));
        cdel = cdel + 1;
    end
    %save(fullfile(dataDir,[expName,'_analysis.mat']),'Conditions','-append')
end


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


%% Ordering PSTH

orderedStr = 'ID ordered';
dans = questdlg('Do you want to order the PSTH other than by IDs?',...
    'Order', 'Yes', 'No', 'No');
ordSubs = 1:nnz(filterIdx(2:Ncl+1));
pclID = gclID(filterIdx(2:end));
if strcmp(dans, 'Yes')
    %     if ~exist('clInfo','var')
    %         clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
    %     end
    % varClass = varfun(@class,clInfo,'OutputFormat','cell');
    [ordSel, iOk] = listdlg('ListString', Data.Properties.VariableNames,...
        'SelectionMode', 'multiple');
    orderedStr = [];
    ordVar = Data.Properties.VariableNames(ordSel);
    for cvar = 1:numel(ordVar)
        orderedStr = [orderedStr, sprintf('%s ',ordVar{cvar})]; %#ok<AGROW>
    end
    orderedStr = [orderedStr, 'ordered'];

    if ~strcmp(ordVar,'id')
        [~,ordSubs] = sortrows(Data(pclID,:),ordVar, 'ascend');
    end
end
%% Plot PSTH
% goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
while size(csNames,1) > size(cst)
    csNames(end) = [];
end

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
        if abs(log10(var(stims(cs,:),[],2))) < 15
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


%% Getting cluster info and adding variables to table

ActiveUnit = false(height(Data),1);
Data = addvars(Data,ActiveUnit,'NewVariableNames','ActiveUnit','After','id');
Data{goods, 'ActiveUnit'} = true;

window = diff(responseWindow);

for a = 1: length(consCondNames)
    Data{Data.ActiveUnit == true,[consCondNames{1,a}, '_Rate_Spont']} = mean(Counts{a,1}')'/window;
    Data{Data.ActiveUnit == true,[consCondNames{1,a}, '_Rate_Evoked']} = mean(Counts{a,2}')'/window;
end


% Significant mechanical responses per condition
b = length(consCondNames);
for a = 1 : length(consCondNames)
    Data{Data.ActiveUnit == true,[consCondNames{1,a}, '_R']} = Results(b).Activity(1).Pvalues < 0.05;
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
            Data{Data.ActiveUnit == true,[consCondNames{1,a},'_vs_',consCondNames{1,b}, '_', Results(e).Activity(c).Type, '_Response']} = Results(e).Activity(c).Pvalues < 0.05;
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
writeClusterInfo(Data, fullfile(figureDir,'cluster_info_TonicResponses.tsv'));

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

med = 0.2;
low = 0.4;
high = 0;
colours = ones(1,3);
colours(1:4,:) = med; colours(5:8,:) = low; colours(9:10,:) = high;

for a = 1%:length(pwrs)
    pwr = pwrs(a);
    % MchTblInd = ['Mech_Control_', num2str(pwr), 'mW_MR'];
    % LasTblInd =  [20,21,22];%['Laser_Control_', num2str(pwr), 'mW_LR'];
    % MchCondControl = ['Mech_Control_', num2str(pwr), 'mW'];
    % LasCondControl = ['Laser_Control_', num2str(pwr), 'mW'];
    % MchLasCond = ['Mech_Laser_', num2str(pwr), 'mW'];
    % EffectTblInd = ['Mech_Control_', num2str(pwr), 'mW_vs_Mech_Laser_', num2str(pwr), 'mW_Evoked_Response'];
    % TblInd = find(clInfo.ActiveUnit); % ATM this only makes rasters that show sig control mech response
    % clIDind = clInfo.id(TblInd);
    %     pwrInd = Power == pwr;
    clIDind =  Data.id(Data.Mech_Low_R);
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
        rasCondSel = [1 2 3];
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
        ax = gobjects(4*Nrcond*Nrcl,1);
        lidx = 1;
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

                stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
                %                 stims = stims([2,3],:);

                stims = stims - median(stims,2);



                for cs = 2 %for cs = 1:size(stims,1)
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




                for cs = 2 %for cs = 1:size(stims,1)
                    if stims(1,cs) > 0.5
                        stim = ones(size(stims(:,cs)))-stims(:,cs);
                    else
                        stim = stims(:,cs);
                    end
                    %                     dP=[0; diff(stim)];
                    stim = smooth(stim,10);
                    ax(lidx) = subplot(4*Nrcond, Nrcl, lidx);
                    if exist('IDs','var')
                        plot(trigTX,stim, 'LineStyle','-','LineWidth', 1,...
                            'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
                        %                         P = stim'; %ax(lidx).Children.YData;
                        %                         dP = [0 diff(P)];
                        %                         dP = smooth(dP,100);
                        %                         dP = dP-mean(dP(1:100));dP=dP/max(dP);
                        %
                        %                         yyaxis right
                        %                         plot(trigTX, dP, 'LineStyle','-','LineWidth', 1,...
                        %                             'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
                    else
                        plot(trigTX,stim,'LineStyle','-','LineWidth',1, 'Color', stmClr(cs,:))

                        %                         P = ax(lidx).Children.YData;
                        %                         dP = [0 diff(P)];
                        %                         dP = smooth(dP,100);
                        %                         dP = dP-mean(dP(1:100));dP=dP/max(dP);
                        %                         yyaxis right
                        %                         plot(trigTX, dP, 'LineStyle','-','LineWidth', 1,...
                        %                             'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
                    end
                    ax(lidx).Visible = 'off';


                    %                                        ax2.Children(1).Color = defineColorForStimuli(IDs(cs));

                    if cs == 1
                        ax2.NextPlot = 'add';
                    end

                end
                %                 ax = gca;
                %
                %
                % %                 ax.YAxis(2).Limits = [0.015, 1];
                % %                 ax.YAxis(2).Visible = 'off';
                %                 ax.FontName ='Arial';
                %                 ax.FontSize = 12;

                %                 f=get(gca,'Children');
                %                 legend(f)
                %
                lidx = lidx + 1;


                %                 lidx = ccl + (cc - 1) * Nrcl;
                ax(lidx) = subplot(4*Nrcond, Nrcl, lidx:lidx+2);       %  subplot( Nrcl, Nrcond, lidx); % to plot the other way around
                title(ax(lidx),sprintf(rasCondNames{cc}), 'Interpreter', 'none') % ,pclID{clSel(ccl)}
                plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                    timeLapse, fs,'',ax(lidx));
                %                    plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                %                     timeLapse, fs,'',colours(lidx,:),ax(lidx));
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

                %                 ax(lidx).XAxis.Visible = 'off';
                ax(lidx).YAxis.Visible = 'off';
                stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
                %                 stims = stims([2,3],:);

                stims = stims - median(stims,2);



                for cs = 2 %for cs = 1:size(stims,1)
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



                lidx = lidx + 3;


            end
        end





        rasConds = rasCondNames{1};
        if length(rasCondNames) > 1
            for r = 2:length(rasCondNames)
                rasConds = [rasConds, '+', rasCondNames{r}];
            end
        end


        ax(2).Title.Color = colours(2,:);
        ax(6).Title.Color = colours(6,:);
        ax(10).Title.Color = colours(10,:);
        ax(1).Children.Color = colours(1,:);
        ax(5).Children.Color = colours(5,:);
        ax(9).Children.Color = colours(9,:);


        linkaxes(ax,'x')
        rasFigName = ['Unit_', cell2mat(cl), '_', ];
        rasFig.Name = [rasFigName, '_', num2str(pwr), 'mW'];
        configureFigureToPDF (rasFig);
        set(rasFig, 'Position', get(0, 'ScreenSize')/2);
        saveas(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds,'_', num2str(timeLapse(1)), '_to_', num2str(timeLapse(2)),'.emf']));
        %savefig(rasFig,fullfile(rasterDir, [rasFigName, ' ', num2str(pwr), 'mW.fig']));
        savefig(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds, '.fig']));
    end
end


%% print mech responses to pdf
targetdir = [figureDir filesep 'difMech']; %put all individual pdfs here, replace blocking with whatever makes sense
mergedir=print_all_figs(targetdir,'-dpdf');  %print all open figures to pdf
cd(figureDir)
merge_PDF_dir(mergedir)  %merge all single figure pdfs to one large pdf, in this directory


%% Ordering PSTH
filterIdx = [true; ismember(gclID, Data.id(Data.Mech_Low_R))];
orderedStr = 'ID ordered';
ordSubs = 1:nnz(filterIdx(2:Ncl+1));
pclID = gclID(filterIdx(2:end));




%% PSTHs for comparisons
binSz = 0.5;
% goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
while size(csNames,1) > size(cst)
    csNames(end) = [];
end

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
    PSTH = PSTH./binSz/sum(delayFlags(:,1));

end

%%
for unit = 1:size(PSTH, 1)
    figure('Name',['unit_', pclID{unit}], 'Color', 'white')
    hold on
    for ccond = 1:length(consideredConditions)
        plot(PSTH(unit,:,ccond));
    end
    ax = gca;
    legend(consCondNames);
    ax.XTickLabel = ax.XTick*binSz-2;

end


%% PopPSTH By Group
figure('Name','DifMech_PopPSTH', 'Color','white'); hold on
[Ncl, Npt, Nconds] = size(PSTH);
psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
for ccond = 1:length(consideredConditions)


    medPSTH = median(sum(PSTH(:,1:40,ccond),1,'omitnan')/(Ncl * sweeps * binSz));
    popPSTH = sum(PSTH(:,:,ccond),1,'omitnan')/(Ncl * sweeps * binSz);
    %         popPSTH = popPSTH-medPSTH;
    popPSTH = smooth(popPSTH, 5);
    plot(psthTX,popPSTH, 'LineWidth',1.5)


end
leg = legend;
leg.String = consCondNames;
leg.Box = 'off';
leg.FontName = 'Arial';
ax = gca;
ax.FontSize = 25;
leg.Location = 'northeast';


%% Absolute Pressure Difference Comparison

cs = 2;
stims = mean(cst,3);
[m,b] = lineariz(stims(cs,:),1,0);
stims(cs,:) = m*stims(cs,:) + b;
% figure('Color','white', 'Name', 'AvPressure');
% plot(stims(cs,:))

figure('Color','white', 'Name', 'Pressures');

hold on

for ccond = 1:Nccond
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    if abs(log10(var(stims(cs,:),[],2))) < 15
        stims(cs,:) = m*stims(cs,:) + b;
        stims(cs,:) = stims(cs,:) - min(stims(cs,:));

        stims(cs,:) = smooth(stims(cs,:),10^4);
    else
        stims(cs,:) = zeros(1,Nt);
    end
    plot(stims(cs,:))
end


%% Delta Pressure Difference Comparison

colours = [0,0,0.75; 0, 0.75, 0; 0.75, 0, 0.75];
cs = 2;
stims = mean(cst,3);
[m,b] = lineariz(stims(cs,:),1,0);
stims(cs,:) = m*stims(cs,:) + b;
% figure('Color','white', 'Name', 'AvPressure');
% plot(stims(cs,:))

figure('Color','white', 'Name', '\Delta Pressures');

hold on

for ccond = 1:Nccond
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    if abs(log10(var(stims(cs,:),[],2))) < 15
        stims(cs,:) = m*stims(cs,:) + b;
        stims(cs,:) = stims(cs,:) - min(stims(cs,:));

        stims(cs,:) = smooth(stims(cs,:),10^4);
    else
        stims(cs,:) = zeros(1,Nt);
    end
    stim = diff(stims(cs,:));
    stim = smooth(stim, 10^4);
    plot(stim, 'color', colours(ccond,:), 'LineWidth', 1.25);
end
fig = gcf;
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 25;

%% Rasters

csNames = fieldnames(Triggers);
% csNames = csNames(2:end);
IDs = csNames;
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));


rasterDir = fullfile(figureDir,'Rasters\');
if ~mkdir(rasterDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end

colours = [0,0,0.75; 0, 0.75, 0; 0.75, 0, 0.75];
cs = 2;

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
    rasCondSel = [1 2 3];
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
    ax = gobjects(4*Nrcond*Nrcl,1);
    lidx = 1;
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

           stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
    if abs(log10(var(stims(cs,:),[],2))) < 15
        stims(cs,:) = m*stims(cs,:) + b;
        stims(cs,:) = stims(cs,:) - min(stims(cs,:));

        stims(cs,:) = smooth(stims(cs,:),10^4);
    else
        stims(cs,:) = zeros(1,Nt);
    end
    stim = [0, diff(stims(cs,:))];
    stim = smooth(stim, 10^4);






            ax(lidx) = subplot(4*Nrcond, Nrcl, lidx);
            if exist('IDs','var')
                plot(trigTX,stims(cs,:), 'LineStyle','-','LineWidth', 1,...
                    'DisplayName', IDs{cs}, 'Color',colours(rasCondSel(cc),:));

            else
                plot(trigTX,stim(cs,:),'LineStyle','-','LineWidth', 1, 'color', colours(rasCondSel(cc),:));

            end
            ax(lidx).Visible = 'off';




            if cs == 1
                ax2.NextPlot = 'add';
            end

        end

        lidx = lidx + 1;


        ax(lidx) = subplot(4*Nrcond, Nrcl, lidx:lidx+2);       %  subplot( Nrcl, Nrcond, lidx); % to plot the other way around
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

        ax(lidx).XAxis.Visible = 'off';
        ax(lidx).YAxis.Visible = 'off';
        stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);


        stims = stims - median(stims,2);





        lidx = lidx + 3;


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
set(rasFig, 'Position', get(0, 'ScreenSize')/2);
saveas(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds,'_', num2str(timeLapse(1)), '_to_', num2str(timeLapse(2)),'.emf']));
%savefig(rasFig,fullfile(rasterDir, [rasFigName, ' ', num2str(pwr), 'mW.fig']));
savefig(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds, '.fig']));
end
