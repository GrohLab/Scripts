function fig = plotUnitInfo(ID, sortedData, clInfo, clWaveforms, Conditions, samplingFrequency, Triggers)



clIDs = ID;
fs = samplingFrequency;
clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.AbsDepth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
Dpth = -Dpth;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end

spkSubs = cellfun(@(x) round(x.*fs), sortedData(spkInd,2),...
'UniformOutput', false);

TriggerTimes = Conditions(1).Triggers(1,:); % assumes Condition(1).Triggers is LaserALL

[firingRatePerCluster, deltaTrigTimeSum, sponSpks, sponIsi] =...
    getSpontFireFreq(spkSubs, TriggerTimes, [0, Conditions(1).Triggers(end,end)], fs, 2);


if class(ID) == 'cell'
    figureName = ID{1};
else
    figureName = ID;
end
% name(strfind(name, '_')) = ' ';
figureName = ['Unit Info: ', figureName];
fig = figure('Name', figureName, 'Color', 'White');



% Preparing condIndices (lIndices)
laser = false(length(Conditions),1);
alltriggers = false(length(Conditions),1);
mech = false(length(Conditions),1);
Hz10 = false(length(Conditions),1);
for a = 1:length(Conditions)
    laser(a) = contains(Conditions(a).name, 'laser', 'IgnoreCase', true);
    alltriggers(a) = contains(Conditions(a).name, 'alltriggers', 'IgnoreCase', true);
    mech(a) = contains(Conditions(a).name, 'mech', 'IgnoreCase', true);
    Hz10(a) = contains(Conditions(a).name, '10Hz', 'IgnoreCase', true);
end
lIndices = find(laser & ~mech & ~alltriggers);
Hz10 = Hz10 & alltriggers;


OptInd = find(Hz10,1,'last');
if ~isempty(OptInd)
    TaggedIDs = TagID([3,6], [0.1,3], clInfo, ID, Conditions((OptInd)).Triggers, sortedData, fs);
else
    TaggedIDs = [];
end
if ~isempty(TaggedIDs)
    rasclr = [0, 0, 1];
else
     rasclr = [0, 0, 0];
end


%% Triggered Dose_Response

subplot(2,2,1)
plotTriggeredDoseResponse(ID, sortedData, lIndices, Conditions, fs); % Need to find the CondIndices!!
ax = gca;

ax.Box = 'off';
ax.FontSize = 9;
ax.FontName = 'Arial';
%% Rasters & PSTHs
% Defaulting to Mech Control for Rasters & PSTH


ctrl = false(length(Conditions),1);
mech = false(length(Conditions),1);
for a = 1:length(Conditions)
    ctrl(a) = contains(Conditions(a).name, 'control', 'IgnoreCase', true);
    mech(a) = contains(Conditions(a).name, 'mech', 'IgnoreCase', true);
end
mcIndices = find(ctrl & mech);
if ~isempty(mcIndices)
    CondTriggers = sort(cat(1,Conditions(mcIndices).Triggers), 'Ascend');
    clr = [1,0,0];
    
    figtitle = 'Mechanical Control Response';
    
else
    figtitle = [{'1Hz Responses (power ascending)'}, {'10Hz Responses (power ascending)'}, {'5Sec Pulse Responses (power ascending)'}];
    for i = 1:3
        CondTriggers = Conditions(2).Triggers(i:3:end,:);
        clr = [0,1,1];
        
        
        
        
        % Rasters
        
        %subplot(6,2,4*i-2)
        subplot(3,2,2*i)
        
        
        plotUnitTriggeredRaster(ID, clInfo, sortedData, CondTriggers, fs, Triggers, clr, rasclr)
        
        ax = gca;
        xlabel(ax,'Time_{(ms)}');
        ax.Title.String = figtitle{i};
        ax.FontSize = 8;
        ax.FontName = 'Arial';
        
        
%         % PSTH
%         
%         subplot(6,2,4*i)
%         
%         plotUnitTriggeredPSTH(ID, clInfo, sortedData, CondTriggers, fs, Triggers, clr)
%         
%         ax = gca;
%         ax.FontSize = 8;
%         ax.FontName = 'Arial';
    end
end
%% Waveform

subplot(4,2,7)

inds = find(ismember(clWaveforms(:,1), ID));
if isempty(inds)
    regionID_end = strfind(clWaveforms{1,1}, '_');
    for a = 1:length(clWaveforms)
        clWaveforms{a,1} = clWaveforms{a,1}(regionID_end+1:end);
    end
    inds = find(ismember(clWaveforms(:,1), ID));
end

ind = inds(1);
[bartho, peak2Ind] = getBarthoMeasure((mean(clWaveforms{ind,2},2)), fs);
plot(mean(clWaveforms{ind,2},2), 'Color', 'b');
ax = gca;

sz = size(clWaveforms{ind,2}(:,1));
sz = sz(1);

 trough = min(mean(clWaveforms{ind,2},2));
% troughInd = find(ismember(mean(clWaveforms{ind,2},2), trough));
% peak2 = max(mean(clWaveforms{ind,2}(troughInd:end,:),2));
% % peak2Ind = find(ismember(mean(clWaveforms{ind,2},2), peak2));
% hold on
% plot(peak2Ind*ones(100), linspace(trough, peak2, 100),...
%     'LineStyle', '--', 'Color', 'k');
% plot(linspace(troughInd, peak2Ind,peak2Ind-troughInd), ...
%     trough*ones(peak2Ind-troughInd), 'Color', 'k');
% text(peak2Ind + 1, trough, [num2str(round(bartho(1)*1e3, 2)), '_{ms}'])
% hold off
%                       FIX THIS!


ax.XLim = [1, sz];
ax.YLim = [round(1.1*trough,-2), round(-1.1*trough,-2)]; % round here
ax.YTick = round((ax.YLim(1):100:ax.YLim(2)));
ax.YTickLabel = round(ax.YTick./1e3,2);
ax.XTick = linspace(1,sz,10);
ax.XTickLabel = round(ax.XTick./fs*1e3,1);
ax.XTickLabelRotation = -45;
% title(['Unit ',clWaveforms{ind,1},' Waveform'], 'Interpreter', 'none');
xlabel(ax,'Time_{(ms)}');
ylabel(ax,'Amplitude_{(mV)}');
WvMeasure.ID = figureName;
WvMeasure.PeakTrough = bartho(1);
WvMeasure.HalfWidth = bartho(2);
difs = [ax.XLim(2)-ax.XLim(1); ax.YLim(2)-ax.YLim(1)];
text(ax, ax.XLim(1)+0.1*difs(1), ax.YLim(2)-0.1*difs(2), ['Depth = ', num2str(Dpth),'{\mum}'], 'Interpreter', 'tex');
ax.Box = 'off';
ax.FontSize = 9;
ax.FontName = 'Arial';
%% ISIs & Bursting

subplot(4,2,5)



[Nr, Nc] = size(sponSpks);
if iscell(sponSpks) && all(cellfun(@isnumeric,sponSpks))
    Nu = max(Nr,Nc);
    % Validation for subscripts or time points.
    SubsFlag = any(cellfun(@all,cellfun(@eq,...
        cellfun(@minus,...
        cellfun(@round,sponSpks,...
        'UniformOutput',false),sponSpks,...
        'UniformOutput',false),repmat({0},Nr,Nc),...
        'UniformOutput',false)));
    mintd = 1/fs;
    if ~SubsFlag
        fprintf(1,'Given time points instead of subscripts.')
        fprintf(1,' Considering time in seconds\n')
        fs = 1;
    end
else
    fprintf(1,'Not sure how to manage this entry...\n')
    fprintf(1,'Please provide either a cell array, a numeric vector, or a ')
    fprintf(1,'logical time series\n')
    return
end
%% Inter-spike intervals
FRpu = zeros(Nu,3);
for cu = 1:Nu
    lisi = log10(diff(sponSpks{cu}./fs));
    if any(isinf(lisi)) 
        fprintf(1,'Cluster %d ',cu);
        
        fprintf(1,' has repeated time points.\n')
        fprintf(1,'These time points will not be considered.\n')
        lisi(isinf(lisi)) = [];
    end
    if isempty(lisi)
        fprintf(1,'Cluster %d ',cu);
        if IDflag
            fprintf(1,'(%s)',ID{cu})
        end
        fprintf(1,' has one time point.\n')
        fprintf(1,'Skipping...\n')
        continue
    end
    FRpu(cu,:) = [1/(10^mean(lisi,'omitnan')),...
        1/(10^median(lisi,'omitnan')),...
        1/(10^mode(lisi))];
    hisi = histogram(lisi,linspace(log10(mintd),3,100));
    cts = hisi.BinCounts;
    bns = (hisi.BinEdges(1:end-1) + hisi.BinEdges(2:end))/2;
    plot(bns,cts./sum(cts),'LineWidth',1);
    ax = gca;
    ax.FontSize = 9;
    ax.XLim = [-4, 2];
    ax.XTick = [-4:2];
    ax.XTickLabel = 10.^ax.XTick * 1e3;
    ax.XTickLabelRotation = -45;
    xlabel(ax,'Time_{(ms)}'); ylabel(ax,'ISI Probability');
    grid(ax,'on')
    Ncts = cts/sum(cts);
    yyaxis('right');plot(bns,cumsum(Ncts),'LineStyle','-', 'Color', [1,0,0])
    ylabel('Cumulative fraction');ax = gca;
    ax.YAxis(2).Limits = [0, 1];
    ax.YAxis(2).Color = [0.1,0.1,0.1];
end

[burstSpkFreq, burstCF, nBursts, nSpikes, eventRatio] = getBurstingMeasures(ID, 4, 100, 2, sortedData, fs);

difs = [ax.XLim(2)-ax.XLim(1); ax.YLim(2)-ax.YLim(1)];
text(ax, ax.XLim(2)-0.3*difs(1), ax.YLim(2)-0.1*difs(2), ['Burst/Spike Ratio = ', num2str(round(eventRatio, 2))]);
ax.Box = 'off';
ax.FontSize = 9;
ax.FontName = 'Arial';
end
