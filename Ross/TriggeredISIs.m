% 26.5.2021




% function = TriggeredISIs(clInfo,sortedData,Conditions,responseWindow,fs,ISIspar,onOffStr)


%% Getting the ActiveUnit ISIs from sortedData
ind = clInfo.ActiveUnit;
gclID = clInfo.id(ind == true);
%ind = ismember(sortedData(:,1), gclID(wruIdx,1));
ind = ismember(sortedData(:,1), gclID);
Ncl = sum(ind);
spkSubs = cellfun(@(x) round(x.*fs), sortedData(ind,2),...
    'UniformOutput', false);
ISIVals = cellfun(@(x) [x(1)/fs; diff(x)/fs], spkSubs,...
    'UniformOutput', 0);
NnzvPcl = cellfun(@numel,ISIVals);
Nnzv = sum(NnzvPcl);
rows = cell2mat(arrayfun(@(x,y) repmat(x,y,1), (1:Ncl)', NnzvPcl,...
    'UniformOutput', 0));
cols = cell2mat(spkSubs);
vals = cell2mat(ISIVals);
ISIspar = sparse(rows, cols, vals);

%% Creating the TrigISIs Struct


ConsConds = Conditions(18:20);
nCond = length(ConsConds);
% sortedData = sortedData(:,1);
for chCond = 1:nCond
    TrigISIs(chCond).name = ConsConds(chCond).name;
    TrigISIs(chCond).Vals(1).name = 'Spontaneous';
    TrigISIs(chCond).Vals(2).name = 'Evoked';
end

%% Adding ISIs to TrigISIs
spontWindow = -flip(respWindow);

for chCond = 1:nCond
    
    
    % contains will give multiple units when looking for e.g. cl45
%     spkLog = StepWaveform.subs2idx(round(sortedData{goods(1),2}*fs),Ns);
    % spkSubs replaces round(sortedData{goods(1),2}*fs) for the rest of the
    % clusters
    % Subscript column vectors for the rest good clusters
    for wIndex = 1:2
        if wIndex == 1
            Window = spontWindow;
        else
            Window = respWindow;
        end
        [~, isiStack] = getStacks(false,ConsConds(chCond).Triggers, onOffStr,...
            Window,fs,fs,[],ISIspar);
        lInda = isiStack > 0; 
        % timelapse becomes spontaneousWindow for pre-trigger, and responseWindow
        % for post
        TrigISIs(chCond).Vals(wIndex).TriggeredIsI = isiStack;
        for histInd = 1: Ncl
            figure('Visible','off');
            hisi = histogram(log10(isiStack(histInd,:,:)), 'BinEdges', log10(0.001):0.01:log10(10));
            TrigISIs(chCond).Vals(wIndex).cts{histInd} = hisi.BinCounts;
            TrigISIs(chCond).Vals(wIndex).bns{histInd} = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
            
            close gcf;
        end
    end
end
%% ISIs and CumISIs
for chCond = 1:nCond
    for a = 1:length(TrigISIs(chCond).Vals(wIndex).cts)
        TrigISIs(chCond).Vals(1).ISI{a} = TrigISIs(chCond).Vals(1).cts{a}./sum(TrigISIs(chCond).Vals(1).cts{a});
        TrigISIs(chCond).Vals(2).ISI{a} = TrigISIs(chCond).Vals(2).cts{a}./sum(TrigISIs(chCond).Vals(2).cts{a});
        TrigISIs(chCond).Vals(1).CumISI{a} = cumsum(TrigISIs(chCond).Vals(1).ISI{a});
        TrigISIs(chCond).Vals(2).CumISI{a} = cumsum(TrigISIs(chCond).Vals(2).ISI{a});
    end
end
%% Saving ISIhist
% save(fullfile(dataDir,[expName,'_', num2str(responseWindow), '_TriggeredISIshistBase10.mat']), 'ISI', 'ConsConds', '-v7.3');