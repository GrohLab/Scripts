% function ISIbox = TriggeredISIs(clInfo,sortedData,Conditions,responseWindow,fs,ISIspar,onOffStr)
% spontaneousWindow = -flip(responseWindow);
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