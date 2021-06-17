    function fig = plotLatencyJitterDepth(clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)
fs = samplingFrequency;
clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.abs_depth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end
name = ConditionName;
name(strfind(name, '_')) = ' ';
name = ['Latency vs Jitter vs Depth: ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
jitter = randi(20,size(Dpth));
Dpth = Dpth + jitter; % adding small random jitter for visualisation
fig = figure('Name', name, 'Color', 'White');

scatter3(mn,sd,Dpth, 50, [0.5, 0.5, 0.5], '*')

pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
zlim([round(round(min(Dpth),-2)/2,-2)*2-200,0]);
zticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 50]);
xticks(0:5:50);
xlabel('Latency_{(ms)}');
ylabel('StdDev_{(ms)}', 'Interpreter', 'tex');
zlabel('Depth_{(\mum)}', 'Interpreter', 'tex');

% ax.XTickLabelRotation = -45;
%ax = gca;
%hold on
%laser = plot([0:pulseWidth], (ax.YLim(1)+25)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;
% ax.FontName = 'Times New Roman';
% Create rectangle
% lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
lgd.String{1} = 'Unit Mean +/- SD';
%lgd.String{2} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

%rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
%'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
end


