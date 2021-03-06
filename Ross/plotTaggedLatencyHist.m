function fig = plotTaggedLatencyHist(latencyCutoffs, sdCutoffs, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)
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
name = ConditionName;
name(strfind(name, '_')) = ' ';
name = ['Latency Histogram: ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');

tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2);

fig = figure('Name', name, 'Color', 'White');

h = histogram(round(mn(~tagged)));
h.FaceColor = [0.5, 0.5, 0.5];
hold on
j = histogram(round(mn(tagged)));
j.FaceColor = [0, 1, 1];
hold off

pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([0,round(mode(round(mn)+5))]);
%yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 50]);
xticks(0:5:50);
xlabel('Latency_{(ms)}');
ylabel('Frequency_{(no. of units)}', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
laserLngth = linspace(0, pulseWidth + 0.5);
hold on
laser = plot(laserLngth, (ax.YLim(2))*ones(length(laserLngth),1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;
% ax.FontName = 'Times New Roman';
% Create rectangle
%lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);


rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth + 0.5, ax.YLim(2) - ax.YLim(1)], ...
'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])

lgd = legend;
lgd.String{1} = ['Units (n=', num2str(sum(~tagged)), ')'];
lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];



end


