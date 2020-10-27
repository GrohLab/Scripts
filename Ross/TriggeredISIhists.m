% function ISIbox = TriggeredISIs(clInfo,sortedData,Conditions,responseWindow,fs,ISIspar,onOffStr)
% spontaneousWindow = -flip(responseWindow);
ConsConds = Conditions(2:end);
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
    elseif contains(ConsConds(ChCond).name, ['ch_Laser', 'Last'], 'IgnoreCase', true)
       name = [ConsConds(ChCond-1).name, '_MR']; % make this more robust to match laser intensity
    end
    
    Ind = find(clInfo.(name));
    ID = clInfo.id(Ind);
    id = find(contains(gclID, ID))';
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
        figure('Visible','off');
        hisi = histogram(log(isiStack), 'BinEdges', log(1/fs):0.01:log(100));
        ISIhist(ChCond).Vals(wIndex).cts = hisi.BinCounts;
        ISIhist(ChCond).Vals(wIndex).bns = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
    end
end
% save(fullfile(dataDir,[expName,'_ISIhist.mat']), 'ISIhist', 'ConsConds', '-v7.3');