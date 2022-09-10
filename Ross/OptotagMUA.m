function [TaggedIDs, fig] = OptotagMUA(figureDir, latencyCutoffs, sdCutoffs, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)



layer = 'L5';
% expNames = {'M5', 'M6', 'M7'};

fs = samplingFrequency;
clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.abs_depth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
depth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end
name = ConditionName;
% name(strfind(name, '_')) = ' ';
name = ['Optotagging_', name];
tm = 2.5e-2;
[Latencies, Fidelities] = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, tm);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
rng('default');
jitter = randi(20,size(depth));
depth = depth + jitter; % adding small jitter for visualisation
depth = -1*depth;

% Removing and highlighting dodgy units from plot
lowestLat = latencyCutoffs(1);
lowestSD = sdCutoffs(1);
dodgyMean = ID(find(mn<lowestLat));
dodgySD = ID(find(sd<lowestSD));
fprintf(['\n The following units have mean latenices lower than ', num2str(lowestLat), 'ms '...
    '\n and have been removed from the plot: \n']);
for a = 1:length(dodgyMean)
    fprintf([num2str(dodgyMean{a}), '\n']);
end
fprintf(['\n The following units have standard devations lower than ', num2str(lowestSD), 'ms '...
    '\n and have been removed from the plot: \n']);
for a = 1:length(dodgySD)
    fprintf([num2str(dodgySD{a}), '\n']);
end

mn(mn<lowestLat) = NaN;
sd(sd<lowestSD) = NaN;
nan = isnan(mn) | isnan(sd);
fd = cell2mat(Fidelities);

latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');

tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2);
%tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2) & depth > -500 & depth < -300;
taggedInd = spkInd(tagged);
TaggedIDs = sortedData(taggedInd,1);
nontagged = ~tagged & ~nan;

minDepth = depth(min(find(tagged)));
maxDepth = depth(max(find(tagged)));

fig = figure('Name', [name, '_pie'], 'Color', 'White');

% subplot(2,2,1)
% subplot(3,2,2)

x = [sum(nontagged), sum(tagged)];
explode = [0, 1];
% labels = {['Tagged Fraction = ',num2str(tg), '%'], ' '};
pie(x, explode)
ax = gca;
ax.FontName = 'Arial';
if sum(tagged) > 0
    ax.Children(2).FaceColor = [0, 0.5, 1];
    ax.Children(4).FaceColor = [0.5 0.5 0.5];
else
    ax.Children(2).FaceColor = [0.5 0.5 0.5];
end
ax.FontSize = 10;
% ax.FontName = 'Times New Roman';
set(fig, 'Position', [1.6810, 0.2794, 0.5, 0.5]*1.0e+03);

% Formatting (text is in weird place if there are no tagged cells)
xAx = ax.XLim(1) - diff(ax.XLim)/6;
yAx = ax.YLim(1) + diff(ax.YLim)/1.75;
if sum(tagged) == 0
    lgd = legend;
    lgd.String{1} = ['Untagged'];
    text(ax,xAx,yAx+1.0,['Selection Criteria:'],'FontSize', 10 );
    text(ax,xAx,yAx+0.75,['Mean Latencies between ' num2str(latencyCutoffs(1)), ' and ',num2str(latencyCutoffs(2)), 'ms'], 'FontSize', 10);
    text(ax,xAx,yAx+0.50,['Standard Deviations between ' num2str(sdCutoffs(1)), ' and ',num2str(sdCutoffs(2)), 'ms'], 'FontSize', 10);
    text(ax,xAx,yAx-0.25,'Opto-tagged Units Displayed:', 'FontSize', 10);
    text(ax,xAx,yAx-0.5,['Depths between ' num2str(minDepth), ' and ',num2str(maxDepth), '\mum'], 'Interpreter', 'tex', 'FontSize', 10);
    lgd.String{2} = ['Opto-tagged'];
    lgd.Location = 'northoutside';
    %lgd.Box = 'off';
    
else
    lgd = legend;
    lgd.String{1} = ['Untagged'];
    text(ax,xAx,yAx+1.25,['Selection Criteria:'],'FontSize', 10 );
    text(ax,xAx,yAx+1,['Mean Latencies between ' num2str(latencyCutoffs(1)), ' and ',num2str(latencyCutoffs(2)), 'ms'], 'FontSize', 10);
    text(ax,xAx,yAx+0.75,['Standard Deviations between ' num2str(sdCutoffs(1)), ' and ',num2str(sdCutoffs(2)), 'ms'], 'FontSize', 10);
    text(ax,xAx,yAx,'Opto-tagged Units Displayed:', 'FontSize', 10);
    text(ax,xAx,yAx-0.25,['Depths between ' num2str(minDepth), ' and ',num2str(maxDepth), '\mum'], 'Interpreter', 'tex', 'FontSize', 10);
    lgd.String{2} = ['Opto-tagged'];
    lgd.Location = 'northoutside';
    %lgd.Box = 'off';
end
set(fig, 'Position', [1.6810, 0.2794, 0.5, 0.5]*1.0e+03);
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.pdf']),'ContentType','vector');
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.emf']),'ContentType','vector');

fig = figure('Name', [name, '_scatter'], 'Color', 'White');

% subplot(2,2,2)
% subplot(3,2,4)

scatter(mn(~tagged),sd(~tagged),12,depth(~tagged),'filled'); %  1+4*fd(~tagged),  scatter(mn(~tagged),sd(~tagged),12, [0.5, 0.5, 0.5], '*')
hold on
scatter(mn(tagged),sd(tagged),12, depth(tagged),'filled'); % 1+4*fd(tagged), 
hold off
pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([0,10]); %ylim([0,round(round(max(sd))+5,-1)]);
%yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 tm*1e3]);
xticks(0:5:tm*1e3);
xlabel('First-Spike Latency [ms]');
ylabel('StdDev [ms]', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
laser = plot([0:pulseWidth], (ax.YLim(2))*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 10;
% ax.FontName = 'Times New Roman';
% Create rectangle
%lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
lgd.String{1} = ['Non-tagged Units (n=', num2str(sum(nontagged)), ')'];
lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
lgd.FontSize = 10;
set(fig, 'Position', [1.6810, 0.2794, 0.5, 0.5]*1.0e+03);
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.pdf']),'ContentType','vector');
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.emf']),'ContentType','vector');


fig = figure('Name', [name, '_depthlatency'], 'Color', 'White');

% subplot(2,2,3)
% subplot(1,2,1)

% errorbar(mn(~tagged),depth(~tagged),sd(~tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', [0.5, 0.5, 0.5], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 0, 'MarkerFaceColor',[0.5,0.5,0.5]);
% hold on
errorbar(mn(tagged),depth(tagged),sd(tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color',[0, 0.5, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 0, 'MarkerFaceColor',[0,0.5,1]);
hold off


pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([-1400, 0]); % ylim([round(round(min(depth),-2)/2,-2)*2-200,0]);
yticks([round(round(min(depth),-2)/2,-2)*2-200:200:0]);
xlim([0 tm*1e3]);
xticks(0:5:tm*1e3);
xlabel('First-Spike Latency [ms]');
ylabel('Unit Depth [\mum]', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
laser = plot([0:pulseWidth], ax.YLim(2)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 10;
% ax.FontName = 'Times New Roman';
% Create rectangle
% lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

lgd = legend;
% lgd.String{1} = (['Unit Mean +/- SD (n=', num2str(sum(nontagged)), ')']);
lgd.String{1} = (['Opto-tagged Units (n=', num2str(sum(tagged)), ')']);
lgd.String{2} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
lgd.FontSize = 10;
set(fig, 'Position', [1.6810, 0.2794, 0.4, 0.5]*1.0e+03);
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.pdf']),'ContentType','vector');
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.emf']),'ContentType','vector');


fig = figure('Name', [name, '_bar'], 'Color', 'White');
% subplot(2,2,4)
% subplot(3,2,6)

h = histogram(round(mn(~tagged)));
h.FaceColor = [0.5, 0.5, 0.5];
hold on
j = histogram(round(mn(tagged)));
j.FaceColor = [0, 0.5, 1];
hold off

pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
% ylim([0,sum(mode(round(mn)+5))]);
%yticks([round(round(min(Dpth),-2)/2,-2)*2-200:200:0]);
xlim([0 tm*1e3]);
xticks(0:5:tm*1e3);
yMax = max([max(h.BinCounts), max(j.BinCounts)]);
ylim([0, round(yMax + 5, -1)]);
xlabel('First-Spike Latency [ms]');
ylabel('Frequency [no. of units]', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
ax.FontName = 'Arial';
laserLngth = linspace(0, pulseWidth + 0.5);
hold on
laser = plot(laserLngth, (ax.YLim(2))*ones(length(laserLngth),1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax.FontSize = 10;
% ax.FontName = 'Times New Roman';
% Create rectangle
%lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);


rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth + 0.5, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])

lgd = legend;
lgd.String{1} = ['Non-tagged Units (n=', num2str(sum(nontagged)), ')'];
lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
lgd.FontSize = 10;

set(fig, 'Position', [1.6810, 0.2794, 0.5, 0.5]*1.0e+03);

idxTagged = ismember(clInfo.id, TaggedIDs);
% if ~any(ismember(clInfo.Properties.VariableNames,'Tagged'))
%     clInfo = addvars(clInfo,idxTagged,'After','ActiveUnit',...
%         'NewVariableNames','Tagged');
% end


exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.pdf']),'ContentType','vector');
exportgraphics(fig, fullfile(figureDir, [layer, '_', fig.Name, '.emf']),'ContentType','vector');
end


