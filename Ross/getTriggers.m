% function [Conditions, Triggers] = getTriggers(expName, dataDir, fs)

%% Assumes that you do your mech (with/without laser) trials in one batch of recordings, ie. not laser (L), mech (M), L, M
    % e.g. LLLLLMMM is fine, as is MMMMMMMMMLL ...this script will
    % currently not work for LTP-style recordings
    
    % Another assumption atm is that a continuous laser pulse > 100ms occurs in every
    % recording in which laser TTL trggers are present, and that the time
    % between trials cannot be less than this pulse length
    
    
    % One more thing - We're assuming mechanical control is always the
    % first condition of a recording with MechTTL pulses
    
% Make sure that any mechanical file you want to be considered has minMech
% nummber of MechTTL pulses
    
    
    
minMech = 10;
%% Finding and rearranging smrx files
iOk = -1;
impStr = 'Rhd';

% Get fs from dataDir - later




smrxFiles = dir(fullfile(dataDir,'*.smrx'));
if isempty(smrxFiles)
    smrxFiles = dir(fullfile(dataDir,'*\*.smrx'));
elseif isempty(smrxFiles)
    fprintf(1,'The given directory contains no .smrx files. Please ')
    fprintf(1,'try again with another folder which do contain .smrx files.\n')
    return
end

cellSmrxFiles = {smrxFiles.name};
Nf = size(smrxFiles,1);

% Selecting files
[incFiles, iok] = listdlg('ListString',cellSmrxFiles(1,:),...
    'SelectionMode','multiple',...
    'PromptString','Select the files to join:',...
    'InitialValue',1:Nf);
if iok
    smrxFiles = smrxFiles(incFiles);
    cellSmrxFiles = {smrxFiles.name};
    Nf = size(smrxFiles,1);
else
    fprintf(1,'Cancelling...\n')
    return
end

% Merging order
fileOrder = (1:Nf)';
defInput = num2cell(num2str(fileOrder));
answr = inputdlg(cellSmrxFiles,'File order',[1, 60],defInput);
nFileOrder = str2double(answr);
nSmrxFiles = smrxFiles;
if ~isempty(answr) && sum(abs(fileOrder - nFileOrder)) ~= 0
    fprintf(1,'Changing file order...\n')
    nSmrxFiles(nFileOrder) = smrxFiles(fileOrder);
    smrxFiles = nSmrxFiles;
else
    fprintf('File order not altered\n')
end
clearvars nSmrxFiles nFileOrder
%% getting ConditionSignals
for i = 1:length(str2num(cell2mat(answr)))
    getConditionSignalsBF(fopen([dataDir, '\', smrxFiles(i).name]));
    CondSigs(i).name = smrxFiles(i).name(1:end-5);
    CondSigs(i).Sig = load([dataDir, '\', CondSigs(i).name, '_CondSig.mat']);
end
%% getting Trigs



% Getting header names and numbers



for i = 1:length(CondSigs)
    fnames{i} = fieldnames(CondSigs(i).Sig(1));
    headerNames{i} = (fnames{i}(contains(fnames{i}, 'head')));
    strHeaders = cell2mat(headerNames{i}); % This will fail if header > 99;
    headerNumbers{i} = strHeaders(:,5:end);
end


nTriggers = [];
for sZ = 1:length(headerNames)
    nTriggers = [nTriggers; length(headerNames{sZ})];
end

mostTriggers = max(nTriggers);
ind = find(nTriggers==mostTriggers, 1);


% Making the Trigs names

for a = 1:mostTriggers
        
        Trigs(a).name = CondSigs(ind).Sig.(headerNames{ind}{a}).title;
        
        if contains(Trigs(a).name, 'Laser') 
            laserInd = a;
        elseif contains(Trigs(a).name, 'MechTTL')
            mchInd = a;
        end
        
    end

% Getting the different laser intensities fromt the CondSigs struct
pwrMissing = false;
for conscondSig = 1:length(CondSigs)
    
     mWfind = strfind(CondSigs(conscondSig).name, 'mW');
    Power{conscondSig} = CondSigs(conscondSig).name(mWfind-2:mWfind+1);
    if isempty(Power{conscondSig})
        pwrMissing = true;
    elseif contains(Power{conscondSig}, '.')
        Power{conscondSig} = CondSigs(conscondSig).name(mWfind-3:mWfind+1);
    elseif contains(Power{conscondSig}, '_') || contains(Power{conscondSig}, ' ')
        Power{conscondSig} = Power{conscondSig}(2:end);
    end
end
if pwrMissing
    fprintf('Could not find all laser intensities - consider renaming condsig files')
end
    


% Adding traces (info), offset, and triggers to Trigs struct

for conscondSig = 1:length(CondSigs)
    
   
    
    nConds = length(headerNames{conscondSig});
    for cc = 1:nConds
        headID = headerNames{conscondSig}{cc};
        mismatchCheck = ismember(Trigs(cc).name,CondSigs(conscondSig).Sig.(headID).title);
        if sum(mismatchCheck) ~= length(mismatchCheck)
            trigType = 1;
            while sum(mismatchCheck) ~= length(mismatchCheck)
                trigType = trigType + 1;
                mismatchCheck = ismember(Trigs(trigType).name, CondSigs(conscondSig).Sig.(headID).title);
            end
        else
            trigType = cc;
        end
        
        
        chanID = ['chan', headerNumbers{conscondSig}(cc,:)];
        Trigs(trigType).info{conscondSig} = CondSigs(conscondSig).Sig.(chanID);
        
        if contains(Trigs(trigType).name, 'MechTTL') || contains(Trigs(trigType).name, 'Laser')
                Trigs(trigType).offset(conscondSig) = length(CondSigs(conscondSig).Sig.(chanID));
                Obj = StepWaveform(Trigs(trigType).info{conscondSig}, fs);
                Trigs(trigType).Triggers{conscondSig} = Obj.subTriggers;
            end
            
            
            
            
            
            
    end
end

%% Starting the Conditions Struct

%% Splitting conditions by laser frequency
cc = 3;
offset = 0;
condHasFreq = [];
trigHasFreq = [];
condHasMech =[];
stimFrequencies = cell(length(Trigs(laserInd).Triggers),1);
pulseLength = cell(length(Trigs(laserInd).Triggers),1);
for a = 1:length(Trigs(laserInd).Triggers)
    if ~isempty(Trigs(mchInd).Triggers{1,a})
                condHasMech = [condHasMech; cc];
    end
    trig = Trigs(laserInd).Triggers{1,a};
    if isempty(trig)
        offset = offset + Trigs(laserInd).offset(a);
    else
        pulseInds = diff(trig') >= 0.1*fs;
        pulsedTriggers = trig(pulseInds,:);
        if ~isempty(pulsedTriggers)
            if range(diff(pulsedTriggers')) >= 0.1*fs
                fprintf('Have you got continuous pulses of different lengths?.\n')
                fprintf('Time to alter the script to sort these.\n')
            end
            pulseLength{a} = round(mean(diff(pulsedTriggers')/fs));
            gaps = sort(unique(round(diff(trig(:,1)),-2))/fs,'ascend');
            minIntervalInd = find(abs(diff([pulseLength{a}*ones(length(gaps),1), gaps]')) < 0.5);
            minIntervalInd = minIntervalInd(1);
            minInt = gaps(minIntervalInd);
            
            stimFrequencies{a} = sort(round(1./gaps(1:minIntervalInd-1)), 'ascend');
        else
            %gaps = sort(unique(round(diff(trig(:,1)),-3))/fs,'ascend'); % What if a condition has no pulses? e.g. Opto-tagging
            gaps = sort(unique(diff(trig(:,1)))/fs,'ascend');
            stimFrequencies{a} = round(1/gaps(1)); % assumption here that condition without continuous pulse has only one condition!
            
        end
        
        for i = 1:length(stimFrequencies{a})
            stimInd = [diff(trig(:,1))> 0.9*fs/stimFrequencies{a}(i) & diff(trig(:,1)) < 1.1*fs/stimFrequencies{a}(i); false];
            shiftInd = [false; stimInd(1:end-1)];
            combInd = stimInd | shiftInd; % come back to this, if last condition to be ran is not continuous pulse , then it might shave off the final trigger
            Conditions(cc).name = ['Laser_', num2str(stimFrequencies{a}(i)), 'Hz_' num2str(Power{a}), '_AllTriggers'];
            Conditions(cc).Triggers = trig(combInd,:) + offset;
            condHasFreq = [condHasFreq, cc];
            trigHasFreq = [trigHasFreq, a];
            cc = cc + 1;
        end
        if ~isempty(pulsedTriggers)
            Conditions(cc).name = ['Laser_', num2str(pulseLength{a}), 'sec_pulse_' num2str(Power{a}), '_AllTriggers']; % Change how this is calculated to better assess trigger allocation[
            Conditions(cc).Triggers = pulsedTriggers + offset;
            minInt = pulseLength{a}*1.1;
            cc = cc + 1;
        end
        
        offset = offset + Trigs(laserInd).offset(a);
    end
    
end
Conditions(1).name = 'Laser_AllTriggers';
Conditions(1).Triggers = sort(cat(1, Conditions.Triggers),'ascend');

% Neatening up Conditions
emptycells = [];
for e = 1:length(Conditions)
    if isempty(Conditions(e).name)
        emptycells = [emptycells; e];
    end
end
emptycells = sort(emptycells, 'descend');
for e = 1:length(emptycells)-1
    emptyInd = emptycells(e);
 Conditions(emptyInd) = [];
end
cc = length(Conditions)+1;
%% Getting laser conditions by block
if ~isempty(condHasFreq)
    blockStart = cc;
    lngth = length(Conditions);
    for i = 3: lngth % can make this condHasFreq if 5sec pulses ae added to AllTriggers_Block
        trigStart = [true; diff(Conditions(i).Triggers(:,1)) >= minInt*fs];
        trigFin = [diff(Conditions(i).Triggers(:,1)) >= minInt*fs; true]; % look at this
        Conditions(cc).name = [Conditions(i).name(1:end-11), 'Block'];
        Conditions(cc).Triggers = [Conditions(i).Triggers(trigStart,1), Conditions(i).Triggers(trigFin,2)];
        if ismember(i, condHasMech)
            condHasMech = [condHasMech; cc];
        end
        cc = cc + 1;
    end
    Conditions(2).name = 'Laser_AllTriggers_Block';
    Conditions(2).Triggers = sort(cat(1, Conditions(blockStart:end).Triggers),'ascend');
else
    Conditions(2) = [];
    cc = cc - 1;
    blockStart = 2;
    condHasMech = condHasMech - 1;
end

%% Splitting mechanical conditions
offset = 0;
mechStart = cc;
trigHasFreq = unique(trigHasFreq);
c = 0;
for a = 1:length(Trigs(mchInd).Triggers)
    if length(Trigs(mchInd).Triggers{1,a}) < minMech % assuming no smrx recording will purposefully contain fewer than minMech TTL pulses
        offset = offset + Trigs(mchInd).offset(a);
    else
        % Assuming mechanical control is first condition
        
        trig = Trigs(mchInd).Triggers{1,a} + offset;
        offset = offset + Trigs(mchInd).offset(a);
        firstLaser = [];
        
        
        nMechConds = length(stimFrequencies{a}) + length(pulseLength{a}) + 1;
         
        for i = 1:nMechConds-1
            c = c + 1;
            ccond = condHasMech(c);
            
            laserPairing = Conditions(ccond).name;
            HzInd = strfind(laserPairing, 'Hz');
            if HzInd ~= false
                Hz = laserPairing(HzInd-2:HzInd-1);
                if contains(Hz, '_') || contains(Hz, ' ')
                    Hz = [Hz(2:end), 'Hz'];
                else
                    Hz = [Hz, 'Hz'];
                end
            else
                HzInd = strfind(laserPairing, 'sec');
                Hz = laserPairing(HzInd-2:HzInd+2);
                if contains(Hz, '_') || contains(Hz, ' ')
                    Hz = Hz(2:end);
                end
            end
            
            Conditions(cc).name = ['Mech_Laser_',Hz, '_' num2str(Power{a})];
            Conditions(cc).Triggers = trig(1+i:nMechConds:end,:);
            cc = cc + 1;
        end
        
     
        Conditions(cc).name = ['Mech_Control_', num2str(Power{a})];
        Conditions(cc).Triggers = trig(1:nMechConds:end,:);
        cc = cc + 1;
    end
 end
if length(Conditions)>= mechStart
    mechEnd = cc;
    Conditions(cc).name = 'Mech_All';
    Conditions(cc).Triggers = sort(cat(1, Conditions(mechStart:end).Triggers),'ascend');
    cc = cc + 1;
end

%% Getting laser control blocks

% This isn't working rn - because you're comparing blockstart with
% mechstart regardless of what blockstart represents...you dummy

if length(Conditions)>= mechStart
    nComparisons = (mechEnd-mechStart)/nMechConds;
    lInd = condHasMech(condHasMech >= blockStart); % Need to sort out condHasMEch!
    mInd = mechStart;
    for a = 1:nComparisons % change to mechEnd so it doesn't do something weird if LaserFreq is last recording
        las = round((Conditions(lInd(a)).Triggers(:,1)),-6);
        mech = round((Conditions(mInd).Triggers(:,1)),-6);
        controlInd = ~ismembertol(las, mech);
        Conditions(cc).name = ['Laser_Control', Conditions(mInd).name(11:end)];
        Conditions(cc).Triggers(:,1) = Conditions(lInd(a)).Triggers(controlInd,1);
        Conditions(cc).Triggers(:,2) = Conditions(lInd(a)).Triggers(controlInd,2);
        mInd = mInd + nMechConds;
        cc = cc + 1;
    end
end


%% Getting Triggers
for a = 1:mostTriggers
    if contains(Trigs(a).name, ' ')
        spc = strfind(Trigs(a).name, ' ');
        Trigs(a).name(spc) = [];
    end
    Triggers.(Trigs(a).name) = [];
    for i = 1:length(Trigs(a).info)
        if isempty(Trigs(a).info{i})
            offset = Trigs(mchInd).offset(i);
            Triggers.(Trigs(a).name) = [Triggers.(Trigs(a).name); zeros(offset,1)];
        else
            
            Triggers.(Trigs(a).name) = [Triggers.(Trigs(a).name); Trigs(a).info{i}];
        end
    end
end

%% Save
save(fullfile(dataDir, [expName, '_analysis.mat']), 'Conditions', 'Triggers', '-v7.3');

% end
