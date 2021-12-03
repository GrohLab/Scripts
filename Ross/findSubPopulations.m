function [TaggedIDs, fig] = findSubPopulations(nkm, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)

% This function plots latency and jitter of latnecy against depth for each
% unit and groups the units according to the measures by using k means
% clustering. nkm determines the number of groups one wants to cluster in
% to.

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
name(strfind(name, '_')) = ' ';
name = ['SubPopulations ', name];
tm = 5e-2;
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, tm);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
rng('default');
jitter = randi(20,size(depth));
depth = depth + jitter; % adding small jitter for visualisation
depth = -1*depth;
kmat = [mn, sd, depth];

 
% Removing and highlighting dodgy units from plot
lowestLat = 2; 
lowestSD = 0.2;
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


km = kmeans(kmat(~nan,:), nkm, 'Replicates', 50);
TaggedIDs = cell(nkm,1);

colours = [0,0,0; 1,0,0; 0,1,0; 0,0,1; 1,1,0; 1,0,1];

fig = figure('Name', name, 'Color', 'White');


for sub = 1:nkm
    TaggedIDs{sub,1} = sub;
    TaggedIDs{sub,2} = ID(km==sub);
    errorbar(mn(km == sub),depth(km == sub),sd(km == sub), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color', colours(sub,:), 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 0, 'MarkerFaceColor', colours(sub,:));
    hold on
end
hold off


legend
lgd = legend;
for sub = 1:nkm
    lgd.String{sub} = ['Group ', num2str(sub), ' mean latency +/- SD'];
end

pulseWidth = round(median(diff(TriggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([round(round(min(depth),-2)/2,-2)*2-200,0]);
yticks([round(round(min(depth),-2)/2,-2)*2-200:200:0]);
xlim([0 tm*1e3]);
xticks(0:5:tm*1e3);
xlabel('First-Spike Latency [ms]');
ylabel('Unit Depth [\mum]', 'Interpreter', 'tex')

ax = gca;
hold on
laser = plot([0:pulseWidth], ax.YLim(2)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 20;
lgd.String{sub+1} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];

% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
lgd.FontSize = 8;

end

