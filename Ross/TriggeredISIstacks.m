% function ISIbox = TriggeredISIs(clInfo,sortedData,Conditions,responseWindow,fs,ISIspar,onOffStr)
% spontaneousWindow = -flip(responseWindow);







ConsConds = Conditions(3:end);
nCond = length(ConsConds);
for ChCond = 1:nCond
    ISIstack(ChCond).name = ConsConds(ChCond).name;
    ISIstack(ChCond).Vals(1).name = 'Spontaneous';
    ISIstack(ChCond).Vals(2).name = 'Evoked';
end

for ChCond = 1:nCond
    name = [ConsConds(ChCond).name, '_MR'];
    
    if contains(ConsConds(ChCond).name, 'Laser_Con', 'IgnoreCase', true)
        name = [ConsConds(ChCond).name, '_LR'];
    elseif contains(ConsConds(ChCond).name, 'ch_Laser', 'IgnoreCase', true)
       name = [ConsConds(ChCond-1).name, '_MR']; % make this more robust to match laser intensity
    end
    
    Ind = find(clInfo.(name));
    ID = clInfo.id(Ind);
    goods = find(ismember(sortedData(:,1), ID))';
    spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % spkSubs replaces round(sortedData{goods(1),2}*fs) for the rest of the
    % clusters
    % Subscript column vectors for the rest good clusters
    for wIndex = 1:2
        if wIndex == 1
            Window = spontaneousWindow;
        else
            Window = responseWindow;
        end      
        [~, isiStack] = getStacks(spkLog,Conditions(ChCond).Triggers, onOffStr,...
            Window,fs,fs,[],ISIspar);
        % timelapse becomes spontaneousWindow for pre-trigger, and responseWindow
        % for post
        
        ISIstack(ChCond).Vals(wIndex).Stack = isiStack;
        
        
    end
end
save(fullfile(dataDir,[expName,'_ISIstacks.mat']), 'ISIstacks', 'ConsConds', '-v7.3');