% function ISIbox = TriggeredISIs(clInfo,sortedData,Conditions,responseWindow,fs,ISIspar,onOffStr)
% spontaneousWindow = -flip(responseWindow);
ConsConds = Conditions(8:end);
nCond = length(ConsConds);
for ChCond = 1:nCond
    ISIbox(ChCond).name = ConsConds(ChCond).name;
    ISIbox(ChCond).Vals(1).name = 'Spontaneous';
    ISIbox(ChCond).Vals(2).name = 'Evoked';
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
    
    ISIbox(ChCond).Vals(1).cts = cell(size(ID));
    ISIbox(ChCond).Vals(2).cts = cell(size(ID));
    ISIbox(ChCond).Vals(1).bns = cell(size(ID));
    ISIbox(ChCond).Vals(2).bns = cell(size(ID));
    
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
            Window,fs,fs,[],ISIspar(1,:));
        % timelapse becomes spontaneousWindow for pre-trigger, and responseWindow
        % for post
        
        
        
        figure('Visible', 'off', 'Color',[1,1,1]);
        hisi = histogram(log(isiStack), 'BinEdges', log(1/fs):0.01:log(100));
        cts = cell(size(goods));
        bns = cell(size(goods));
        ISIbox(ChCond).Vals(wIndex).cts{1} = hisi.BinCounts;
        ISIbox(ChCond).Vals(wIndex).bns{1} = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
        for a = 2:length(ID)
            
            [~, isiStack] = getStacks(spkLog,Conditions(ChCond).Triggers, onOffStr,...
                Window,fs,fs,[],ISIspar(a,:));
            figure('Visible','off');
            hisi = histogram(log(isiStack), 'BinEdges', log(1/fs):0.01:log(100));
            ISIbox(ChCond).Vals(wIndex).cts{a} = hisi.BinCounts;
            ISIbox(ChCond).Vals(wIndex).bns{a} = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
        end
    end
end
 save(fullfile(dataDir,[expName,'_ISIbox.mat']), 'ISIbox', 'ConsConds', '-v7.3');
        %%
        figure('Color',[1,1,1]);
        mintd = 1/fs;
       % hisi = histogram(log(isiStack),linspace(log10(mintd),3,100), 'DisplayStyle', 'stairs');
        plot(bns{1},cts{1}./sum(cts{1}),'LineWidth',1);
        Ncts = cts{1}/sum(cts{1});
        % yyaxis('right');plot(bns{1},cumsum(Ncts),'LineStyle','-')
        hold on
        for a = 2:length(ID)
            plot(bns{a},cts{a}./sum(cts{a}),'LineWidth',1);
            Ncts = cts{a}/sum(cts{a});
            % yyaxis('right');plot(bns{a},cumsum(Ncts),'LineStyle','-')
        end
        fig = gcf;
        ax = fig.Children;
        ax.XTickLabel = round(exp(cellfun(@str2double,ax.XTickLabel)) * 1e3);
        xlabel(ax,'Time [ms]'); ylabel(ax,'ISI Probability');
        grid(ax,'on')
%         
%         
%         
        
        
        
        % % figure; histogram(log(isiStack(1,:)), 'BinEdges', log(1/fs):0.01:log(100), 'Normalization', 'cumcount')
        %
        %
        %
        %
        % hold on
        % for a = 2:length(goods)
        %     spkLog = StepWaveform.subs2idx(spkSubs{a,:},Ns);
        %
        %     [~, isiStack] = getStacks(spkLog,Conditions(chCond).Triggers, onOffStr,...
        %         timeLapse,fs,fs,[],ISIspar(a,:));
        %     histogram(log(isiStack), 'BinEdges', log(1/fs):0.01:log(100))
        % end
        %
        %
        %
        %
        %
        
        
