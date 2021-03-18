% function [Conditions, Triggers] = getTriggers(expName, dataDir, fs)
%% Finding and rearranging smrx files
iOk = -1;
impStr = 'Rhd';


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
answr =  cellfun(@(x) str2num(x),answr);
for i = 1:length(answr)
    getConditionSignalsBF(fopen([dataDir, '\', smrxFiles(i).name]));
    CondSig(i).name = smrxFiles(i).name(1:end-5);
    CondSig(i).Sig = load([dataDir, '\', CondSig(i).name, '_CondSig.mat']);
end
%% getting Trigs
%  fsName = CondSig(1).name(1:4);
% if contains(fsName, '_') || contains(fsName, ' ')
%         fsName = fsName(1:end-1);
%     end
%
%  fs = load([dataDir, '\', fsName, '_sampling_frequency.mat']);

for i = 1:length(CondSig)
    
    
    fnames = fieldnames(CondSig(i).Sig(1));
    
    nTriggers = sum(contains(fnames, 'head'));
    mWfind = strfind(CondSig(i).name, 'mW');
    Power{i} = CondSig(i).name(mWfind-2:mWfind+1);
    if contains(Power{i}, '_') || contains(Power{i}, ' ')
        Power{i} = Power{i}(2:end);
    end
    headerNames = cell2mat(fnames(contains(fnames, 'head')));
    headerNumbers = headerNames(:,5:end);
    for a = 1:nTriggers
        chanID = ['chan', headerNumbers(a,:)];
        Trigs(a).name = CondSig(i).Sig.(headerNames(a,:)).title;
        Trigs(a).info{i} = CondSig(i).Sig.(chanID);
        if contains(Trigs(a).name, ' ')
            spc = strfind(Trigs(a).name, ' ');
            Trigs(a).name(spc) = [];
            
        elseif contains(Trigs(a).name, 'MechTTL') || contains(Trigs(a).name, 'Laser')
            Trigs(a).offset(i) = length(CondSig(i).Sig.(chanID));
            Obj = StepWaveform(Trigs(a).info{i}, fs);
            Trigs(a).Triggers{i} = Obj.subTriggers;
        end
        if contains(Trigs(a).name, 'Laser')
            laserInd = a;
        elseif contains(Trigs(a).name, 'MechTTL')
            mchInd = a;
        end
    end
end

%% Splitting conditions by laser frequency
cc = 3;
offset = 0;
for a = 1:length(Trigs(laserInd).Triggers)
    
    trig = Trigs(laserInd).Triggers{1,a};
    
    pulseInds = diff(trig') > 0.1*fs;
    pulsedTriggers = trig(pulseInds,:);
    if range(diff(pulsedTriggers')) > 0.1*fs
        fprintf('Have you got continuous pulses of different lengths?.\n')
        fprintf('Time to alter the script to sort these.\n')
    end
    pulseLength = round(mean(diff(pulsedTriggers')/fs));
    gaps = sort(unique(round(diff(trig(:,1)),-3))/fs,'ascend');
    minIntervalInd = find(abs(diff([pulseLength*ones(length(gaps),1), gaps]')) < 0.5);
    
    % assumes that gap between stimuli cannot be shorter than any gap within
    % a stimulus
    stimFrequencies = sort(round(1./gaps(1:minIntervalInd-1)), 'ascend');
    if isnan(pulseLength)
        stimFrequencies = sort(round(1./gaps(1:find(gaps==max(gaps))-1)), 'ascend');
    end
    for i = 1:length(stimFrequencies)
        stimInd = diff(trig(:,1))> 0.9*fs/stimFrequencies(i) & diff(trig(:,1)) < 1.1*fs/stimFrequencies(i);
        shifted = [0; stimInd(1:end-1)];
        stimInd = stimInd | shifted;
         
        Conditions(cc).name = ['Laser_', num2str(stimFrequencies(i)), 'Hz_' num2str(Power{a}), '_AllTriggers'];
        Conditions(cc).Triggers = trig(stimInd,:) + offset;
        cc = cc + 1;
    end
    Conditions(cc).name = ['Laser_', num2str(pulseLength), 'sec_pulse_' num2str(Power{a}), '_AllTriggers'];
    Conditions(cc).Triggers = pulsedTriggers + offset;
    cc = cc + 1;
    offset = offset + Trigs(laserInd).offset(a);
end
Conditions(1).name = 'Laser_AllTriggers';
Conditions(1).Triggers = sort(cat(1, Conditions.Triggers),'ascend');
%% Getting laser conditions by block
if ~isempty(stimFrequencies)
    blockStart = cc;
    intInd = gaps(minIntervalInd);
    for i = 3: length(Conditions)
        trigStart = [true; diff(Conditions(i).Triggers(:,1)) >= intInd*fs];
        trigFin = [diff(Conditions(i).Triggers(:,1)) >= intInd*fs; true];
        Conditions(cc).name = [Conditions(i).name(1:end-11), 'Block'];
        Conditions(cc).Triggers = [Conditions(i).Triggers(trigStart,1), Conditions(i).Triggers(trigFin,2)];
        cc = cc + 1;
    end
    Conditions(2).name = 'Laser_AllTriggers_Block';
    Conditions(2).Triggers = sort(cat(1, Conditions(blockStart:end).Triggers),'ascend');
else
    Conditions(2) = [];
    cc = cc - 1;
    blockStart = 2;
end

% %% Splitting mechanical conditions
% offset = 0;
% mechStart = cc;
% 
% 
% for a = 1:length(Trigs(mchInd).Triggers)
%     if ~isempty(Trigs(mchInd).Triggers{1,a})
%         % Assuming mechanical control is first condition
%         nMechConds = length(stimFrequencies) + length(pulseLength) + 1;
%         trig = Trigs(mchInd).Triggers{1,a} + offset;
%         
%         firstLaser = [];
%         if ~isempty(stimFrequencies)
%             for i = 1:nMechConds-1
%                 firstLaser = [firstLaser; Conditions(2 + i).Triggers(1,1)];
%             end
%         else
%             for i = 1:nMechConds-1
%                 firstLaser = [firstLaser; Conditions(1 + i).Triggers(1,1)];
%             end
%         end
%         
%         sortedFirsts = sort(firstLaser, 'ascend');
%         
%         
%         for i = 1:length(sortedFirsts)
%             if ~isempty(stimFrequencies)
%                 ind = find(ismember(firstLaser, sortedFirsts(i))) + 2;
%             else
%                 ind = find(ismember(firstLaser, sortedFirsts(i))) + 1;
%             end
%             laserPairing = Conditions(ind).name;
%             HzInd = strfind(laserPairing, 'Hz');
%             if HzInd == true
%                 Hz = laserPairing(HzInd-2:HzInd-1);
%                 if contains(Hz, '_') || contains(Hz, ' ')
%                     Hz = [Hz(2:end), 'Hz'];
%                 end
%             else
%                 HzInd = strfind(laserPairing, 'sec');
%                 Hz = laserPairing(HzInd-2:HzInd+2);
%                 if contains(Hz, '_') || contains(Hz, ' ')
%                     Hz = Hz(2:end);
%                 end
%             end
%             
%             PwrInd = strfind(laserPairing, 'mW');
%             Pwr = laserPairing(PwrInd-2:PwrInd-1);
%             if contains(Pwr, '_') || contains(Pwr, ' ')
%                 Pwr = Pwr(2:end);
%             end
%             Conditions(cc).name = ['Mech_Laser_',Hz, '_' num2str(Pwr), 'mW'];
%             Conditions(cc).Triggers = trig(1+i:nMechConds:end,:);
%             cc = cc + 1;
%         end
%         Conditions(cc).name = ['Mech_Control_', num2str(Pwr), 'mW'];
%         Conditions(cc).Triggers = trig(1:nMechConds:end,:);
%         cc = cc + 1;
%     end
%     offset = offset + Trigs(mchInd).offset(a);
% end
% if ~isempty(cat(1, Conditions(mechStart:end).Triggers))
%     Conditions(cc).name = 'Mech_All';
%     Conditions(cc).Triggers = sort(cat(1, Conditions(mechStart:end).Triggers),'ascend');
%     cc = cc + 1;
% end
% %% Getting laser control blocks
% if length(Conditions) >= mechStart
%     nComparisons = mechStart-blockStart;
%     lInd = blockStart;
%     mInd = mechStart;
%     for a = 1:nComparisons
%         las = round((Conditions(lInd).Triggers),-2);
%         mech = round((Conditions(mInd).Triggers),-2);
%         controlInd = ~ismember(las, mech);
%         Conditions(cc).name = ['Laser_Control', Conditions(mInd).name(11:end)];
%         Conditions(cc).Triggers(:,1) = Conditions(lInd).Triggers(controlInd(:,1),1);
%         Conditions(cc).Triggers(:,2) = Conditions(lInd).Triggers(controlInd(:,2),2);
%         lInd = lInd + 1;
%         mInd = mInd + nMechConds;
%         cc = cc + 1;
%     end
% end
% 


%% Getting Triggers
for a = 1:nTriggers
    Triggers.(Trigs(a).name) = [];
    for i = 1:length(Trigs(a).info)
        Triggers.(Trigs(a).name) = [Triggers.(Trigs(a).name); Trigs(a).info{i}];
    end
end

%% Save
% save(fullfile(dataDir, [expName, '_analysis.mat']), 'Conditions', 'Triggers', '-v7.3');

% end
