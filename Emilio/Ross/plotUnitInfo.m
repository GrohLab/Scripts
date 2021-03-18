function fig = plotUnitInfo(ID, sortedData, clInfo, clWaveforms, Conditions, samplingFrequency)



clIDs = ID;
fs = samplingFrequency;
clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.AbsDepth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end

spkSubs = cellfun(@(x) round(x.*fs), sortedData(spkInd,2),...
'UniformOutput', false);

TriggerTimes = Conditions(1).Triggers; % assumes Condition(1).Triggers is LaserALL

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


%% Latency vs Depth


Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, 5e-2);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
jitter = randi(20,size(Dpth));
Dpth = Dpth + jitter; % adding small random jitter for visualisation
subplot(2,2,1)

errorbar(mn,Dpth,sd, 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 5);

pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([round(round(min(Dpth),-2)/2,-2)*2-200,0]);
yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 50]);
xticks(0:5:50);
xlabel('Latency_{(ms)}');
ylabel('Depth_{(\mum)}', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
laser = plot([0:pulseWidth], ax.YLim(2)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontSize = 15;
% ax.FontName = 'Times New Roman';
% Create rectangle
% lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
lgd.String{1} = 'Unit Mean +/- SD';
lgd.String{2} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
lgd.FontSize = 8;

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])

%% Rasters
subplot(2,2,2)




%% Waveform
subplot(2,2,3)

inds = find(ismember(clWaveforms(:,1), ID));
if isempty(inds)
    regionID_end = strfind(clWaveforms{1,1}, '_');
    for a = 1:length(clWaveforms)
        clWaveforms{a,1} = clWaveforms{a,1}(regionID_end+1:end);
    end
    inds = find(ismember(clWaveforms(:,1), ID));
end

ind = inds(1);
bartho = getBarthoMeasure((mean(clWaveforms{ind,2},2)), fs);
plot(mean(clWaveforms{ind,2},2), 'Color', 'b');
ax = gca;
sz = size(clWaveforms{ind,2}(:,1));
sz = sz(1);

trough = min(mean(clWaveforms{ind,2},2));
troughInd = find(ismember(mean(clWaveforms{ind,2},2), trough));
peak2 = max(mean(clWaveforms{ind,2}(troughInd:end,:),2));
peak2Ind = find(ismember(mean(clWaveforms{ind,2},2), peak2));
hold on
plot(peak2Ind*ones(100), linspace(trough, peak2, 100),...
    'LineStyle', '--', 'Color', 'k');
plot(linspace(troughInd, peak2Ind,peak2Ind-troughInd), ...
    trough*ones(peak2Ind-troughInd), 'Color', 'k');
text(peak2Ind + 1, trough, [num2str(round(bartho(1)*1e3, 3)), '_{ms}'])
hold off


ax.FontSize = 15;
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
WvMeasure(a).ID = figureName;
WvMeasure(a).PeakTrough = bartho(1);
WvMeasure(a).HalfWidth = bartho(2);



%% ISIs & Bursting

subplot(2,2,4)



[Nr, Nc] = size(spkSubs);
if iscell(spkSubs) && all(cellfun(@isnumeric,spkSubs))
    Nu = max(Nr,Nc);
    % Validation for subscripts or time points.
    SubsFlag = any(cellfun(@all,cellfun(@eq,...
        cellfun(@minus,...
        cellfun(@round,spkSubs,...
        'UniformOutput',false),spkSubs,...
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
    lisi = log10(diff(spkSubs{cu}./fs));
    if any(isinf(lisi)) 
        fprintf(1,'Cluster %d ',cu);
        if IDflag
            fprintf(1,'(%s)',ID{cu})
        end
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
    ax.FontSize = 15;
    ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;
    xlabel(ax,'Time_{(ms)}'); ylabel(ax,'ISI Probability');
    grid(ax,'on')
    Ncts = cts/sum(cts);
    yyaxis('right');plot(bns,cumsum(Ncts),'LineStyle','-', 'Color', [1,0,0])
    ylabel('Cumulative fraction');ax = gca;
    ax.YAxis(2).Limits = [0, 1];
    ax.YAxis(2).Color = [0.1,0.1,0.1];
end

[burstSpkFreq, burstCF, nBursts, nSpikes, eventRatio] = getBurstingMeasures(ID, 4, 100, 2, sortedData, fs);

text(ax, 0.5, 0.5, ['Burst/SpikeRatio = ', num2str(round(eventRatio, 2))]);

end
