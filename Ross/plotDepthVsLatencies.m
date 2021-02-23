function fig = plotDepthVsLatencies(clInfo, clIDs, TriggerTimes, sortedData, samplingFrequency)

clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.AbsDepth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end

Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, samplingFrequency);
md = (cellfun(@median, Latencies)*1e3);
Dpth = -1*Dpth;
fig = figure('Name', 'Depth vs Latency', 'Color', 'White');

semilogx(md,Dpth, 'LineStyle', 'none', 'Marker', '*', 'Color', 'b');


%text(median(Latencies{i}),Dpth(i),[' ', ID{i}], 'FontSize', 8);
ylim([-1500,0]);
yticks([-1500:100:0]);
xlim([1 50]);
xticks([1:1:10, 15, 20, 30, 40 ,50]);
xlabel('Median Trigger Latency_{(ms)}');
ylabel('Depth_{(um)}')
ax = gca;
% ax.XTickLabelRotation = -45;
ax.FontSize = 20;
ax.FontName = 'Times New Roman';
% Create rectangle
rectangle(ax, 'Position', [1e-7, ax.YLim(1), 10 - 1e-7, ax.YLim(2) - ax.YLim(1)], ...
'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
end


