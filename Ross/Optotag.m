function [fig, TaggedIDs] = Optotag(latencyCutoffs, sdCutoffs, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)
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
name = ['Optotagging: ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, 5e-2);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
nan = isnan(mn);
latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');

tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2);
taggedInd = spkInd(tagged);
TaggedIDs = sortedData(taggedInd,1);
nontagged = ~tagged & ~nan;

minDepth = Dpth(min(find(tagged)));
maxDepth = Dpth(max(find(tagged)));

fig = figure('Name', name, 'Color', 'White');

subplot(2,2,1)
x = [sum(tagged), sum(~tagged)];
explode = [1, 0];
% labels = {['Tagged Fraction = ',num2str(tg), '%'], ' '};
pie(x, explode)
ax = gca;
ax.Children(2).FaceColor = [0.5 0.5 0.5];
ax.Children(4).FaceColor = [0 1 1];
ax.FontSize = 30;
% ax.FontName = 'Times New Roman';

lgd = legend;
lgd.String{1} = ['Opto-tagged'];
text(-4.25,1.0,['Selection Criteria:'],'FontSize', 15 );
text(-4.25,0.75,['Mean Latencies between ' num2str(latencyCutoffs(1)), ' and ',num2str(latencyCutoffs(2)), 'ms'], 'FontSize', 10);
text(-4.25,0.50,['Standard Deviations between ' num2str(sdCutoffs(1)), ' and ',num2str(sdCutoffs(2)), 'ms'], 'FontSize', 10);
text(-4.25,-0.25,'Opto-tagged Units Displayed:', 'FontSize', 15);
text(-4.25,-0.5,['Depths between ' num2str(minDepth), ' and ',num2str(maxDepth), '\mum'], 'Interpreter', 'tex', 'FontSize', 10);
lgd.String{2} = ['Untagged '];
lgd.Location = 'northoutside';
%lgd.Box = 'off';

subplot(2,2,2)

scatter(mn(~tagged),sd(~tagged), 50, [0.5, 0.5, 0.5], '*')
hold on
scatter(mn(tagged),sd(tagged), 50, [0.5, 1, 1], '*')
hold off
pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([0,round(round(max(sd))+5,-1)]);
%yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 50]);
xticks(0:5:50);
xlabel('Latency_{(ms)}');
ylabel('StdDev_{(ms)}', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
laser = plot([0:pulseWidth], (ax.YLim(2))*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontSize = 20;
% ax.FontName = 'Times New Roman';
% Create rectangle
%lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
lgd.String{1} = ['Units (n=', num2str(sum(nontagged)), ')'];
lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
lgd.FontSize = 8;


subplot(2,2,3)

errorbar(mn(~tagged),Dpth(~tagged),sd(~tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 5);
hold on
errorbar(mn(tagged),Dpth(tagged),sd(tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', [0.5, 1, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 5);
hold off


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
ax.FontSize = 20;
% ax.FontName = 'Times New Roman';
% Create rectangle
% lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
lgd.String{1} = (['Unit Mean +/- SD (n=', num2str(sum(nontagged)), ')']);
lgd.String{2} = (['Opto-tagged Units (n=', num2str(sum(tagged)), ')']);
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
lgd.FontSize = 8;


subplot(2,2,4)

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
lgd.String{1} = ['Units (n=', num2str(sum(nontagged)), ')'];
lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
lgd.FontSize = 8;


end


