function fig = plotDepthVsLatencies(clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)

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
% name(strfind(name, '_')) = ' ';
name = ['Depth vs Latency ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, samplingFrequency);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
fig = figure('Name', name, 'Color', 'White');

errorbar(mn,Dpth,sd, 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', 'k', 'MarkerSize', 2.5);




%text(median(Latencies{i}),Dpth(i),[' ', ID{i}], 'FontSize', 8);
title(name, 'Interpreter', 'none');
ylim([round(round(min(Dpth),-2)/2,-2)*2-200,0]);
yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 50]);
xticks(0:5:50);
xlabel('Latency_{(ms)}');
ylabel('Depth_{(um)}')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
plot([0:10], (ax.YLim(1)+25)*ones(11,1), 'Color', [0 1 1 0.5], 'LineWidth', 3)
ax = gca;
ax.FontSize = 20;
% ax.FontName = 'Times New Roman';
% Create rectangle
legend({'mean +/- SD'});

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])
text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);
rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), 10, ax.YLim(2) - ax.YLim(1)], ...
'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
end


