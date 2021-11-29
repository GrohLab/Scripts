%% Separating the populations by depth, latency, jitter, and waveform

name = 'Population Separations';

triggerTimes = Conditions(14).Triggers;
clInd = ismember(clInfo.id, gclID);
depths = table(clInfo.id(clInd), clInfo.abs_depth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end
% name = ConditionName;
% name(strfind(name, '_')) = ' ';
% name = ['Latency vs Depth: ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), triggerTimes, fs, 5e-2);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
rng('default');
jitter = randi(20,size(Dpth));
Dpth = Dpth + jitter; % adding small random jitter for visualisation
nan = isnan(mn);

% latencyCutoffs = sort(latencyCutoffs, 'ascend');
% sdCutoffs = sort(sdCutoffs, 'ascend');

% Sorting by depth, latency, and jitter

tagged = mn > 2 & mn <= 10 & sd >=0.2 & sd <= 4; % Could make these a Nby2 matrix 2 assign groups based on latency, jitter, depth
later = mn > 10 & mn <= 15 & sd >=0.2 & sd <= 4;
later = later & Dpth >-800;
latest = mn > 15 & mn <= 50 & sd >=0.2 & sd <= 50;
latest = latest & Dpth <-700;

other = ~tagged & ~later & ~latest;

nTagged = sum(tagged);
nLater = sum(later);
nLatest = sum(latest);
nOther = sum(other);




%Sorting by waveform (put. FS or WS)

% taggedFS = tagged & fsFlag;
% taggedRS = tagged &wsFlag;
% laterFS = later & fsFlag;
% laterRS = later & wsFlag;
% latestFS = latest & fsFlag;
% latestRS = latest & wsFlag;

%% Plotting

red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];


%Plotting the first data point of each group member for the legend
firstLater = find(later, 1, 'first');
firstLatest = find(latest, 1, 'first');
firstL6 = find(tagged, 1, 'first');

tagged(firstL6) = false;
later(firstLater) = false;
latest(firstLatest) = false;


fig = figure('Name', name, 'Color', 'White');

hold on

errorbar(mn(firstLater),Dpth(firstLater),sd(firstLater), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', red, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', red);
errorbar(mn(firstLatest),Dpth(firstLatest),sd(firstLatest), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', green, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', green);
errorbar(mn(firstL6),Dpth(firstL6),sd(firstL6), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', blue, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', blue);

leg = legend;



hold on

errorbar(mn(later),Dpth(later),sd(later), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', red, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', red);
%errorbar(mn(laterRS),Dpth(laterRS),sd(laterRS), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', [1, 0, 0], 'MarkerSize', 6, 'LineWidth', 0.01, 'CapSize', 5); %, 'MarkerFaceColor', [1, 0, 0]);

errorbar(mn(latest),Dpth(latest),sd(latest), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', green, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', green);
%errorbar(mn(latestRS),Dpth(latestRS),sd(latestRS), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', [0, 1, 0], 'MarkerSize', 6, 'LineWidth', 0.01, 'CapSize', 5); %, 'MarkerFaceColor', [0, 1, 0]);

errorbar(mn(tagged),Dpth(tagged),sd(tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', blue, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', blue);
%errorbar(mn(taggedRS),Dpth(taggedRS),sd(taggedRS), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', [0.25, 0.5, 1], 'MarkerSize', 6, 'LineWidth', 0.01, 'CapSize', 5); %, 'MarkerFaceColor', [0.25, 0.5, 1]);

leg.String(4:end) = [];
leg.String(4:end-1) = [];
leg.Box = 'off';
legend('Location', 'southeast');


pulseWidth = round(median(diff(triggerTimes'/fs)')*1e3);
% title(name, 'Interpreter', 'none');
ylim([round(round(min(Dpth),-2)/2,-2)*2-200,0]);
yticks([-2000:500:0]);
xlim([0 50]);
xticks(0:10:50);
xlabel('Latency [ms]');
ylabel('Depth [\mum]', 'Interpreter', 'tex')
% ax.XTickLabelRotation = -45;
ax = gca;
hold on
laser = plot([0:pulseWidth], ax.YLim(2)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
leg = legend;
leg.String = [{'Group 1 Mean +/- SD'}, {'Group 2 Mean +/- SD'}, {'L6 Mean +/- SD'},  {['Laser Pulse (',num2str(pulseWidth), 'ms)']}];
leg.FontSize = 10;
ax.FontName = 'Arial';
ax.FontSize = 15;
% ax.FontName = 'Times New Roman';
% Create rectangle
% lsText = text(3, ax.YLim(1)+45, 'Laser Pulse', 'FontSize', 10);

% lgd = legend;
% lgd.String{1} = (['Unit Mean +/- SD (n=', num2str(sum(nontagged)), ')']);
% lgd.String{2} = (['Opto-tagged Units (n=', num2str(sum(tagged)), ')']);
% lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
% rectangle(ax, 'Position', [1, ax.YLim(2)-100, 4, 50], ...
% 'LineStyle', 'none', 'FaceColor', [0 1 1 0.5])

rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
    'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])

hold off


%% Pie Chart

x = [nLater, nLatest, nTagged, nOther];
explode = ones(1,4);
labels = [{'Group 1'}, {'Group 2'}, {'L6'}, {'Other'}];

figure('Color', 'White');
pie(x, explode)
leg = legend;
leg.String = labels;
leg.Box = 'Off';
leg.Location = 'eastoutside';
leg.FontSize = 8;


ax = gca;
ax.Children(2).FaceColor = 'yellow';
ax.Children(4).FaceColor = blue;
ax.Children(6).FaceColor = green;
ax.Children(8).FaceColor = red;
ax.FontName = 'Arial';

for slice = 1:2:2*length(x)
ax.Children(slice).FontName = 'Arial';
ax.Children(1).FontSize = 10;
end

ax.Children(2).FaceColor = 'yellow';
ax.Children(4).FaceColor = blue;
ax.Children(6).FaceColor = green;
ax.Children(8).FaceColor = red;

for slice = 2:2:2*length(x)
    ax.Children(slice).EdgeColor = [1,1,1];
end
