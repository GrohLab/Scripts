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

L6 = mn > 2 & mn <= 10 & sd >=0.2 & sd <= 4; % Could make these a Nby2 matrix 2 assign groups based on latency, jitter, depth
L6 = L6 & Dpth <-650;
L4 = mn > 10 & mn <= 15 & sd >=0.5 & sd <= 10;
L4 = L4 & Dpth >-700;
L5 = mn > 15 & mn <= 50 & sd >=0.2 & sd <= 50;
L5 = L5 & Dpth >-900 & Dpth <-600;

other = ~L6 & ~L4 & ~L5;

nL6 = sum(L6);
nL4 = sum(L4);
nL5 = sum(L5);
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
firstL4 = find(L4, 1, 'first');
firstL5 = find(L5, 1, 'first');
firstL6 = find(L6, 1, 'first');

L6(firstL6) = false;
L4(firstL4) = false;
L5(firstL5) = false;


fig = figure('Name', name, 'Color', 'White');

hold on

errorbar(mn(firstL4),Dpth(firstL4),sd(firstL4), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', green, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', green);
errorbar(mn(firstL5),Dpth(firstL5),sd(firstL5), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', red, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', red);
errorbar(mn(firstL6),Dpth(firstL6),sd(firstL6), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', blue, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', blue);

leg = legend;



hold on

errorbar(mn(L4),Dpth(L4),sd(L4), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', green, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', green);
%errorbar(mn(laterRS),Dpth(laterRS),sd(laterRS), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', [1, 0, 0], 'MarkerSize', 6, 'LineWidth', 0.01, 'CapSize', 5); %, 'MarkerFaceColor', [1, 0, 0]);

errorbar(mn(L5),Dpth(L5),sd(L5), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', red, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', red);
%errorbar(mn(latestRS),Dpth(latestRS),sd(latestRS), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', [0, 1, 0], 'MarkerSize', 6, 'LineWidth', 0.01, 'CapSize', 5); %, 'MarkerFaceColor', [0, 1, 0]);

errorbar(mn(L6),Dpth(L6),sd(L6), 'horizontal', 'LineStyle', 'none', 'Marker', 'o', 'Color', blue, 'MarkerSize', 2, 'LineWidth', 0.01, 'CapSize', 5, 'MarkerFaceColor', blue);
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
leg.String = [{'Put L4 Mean +/- SD'}, {'L5 Mean +/- SD'}, {'Put L6 Mean +/- SD'},  {['Laser Pulse (',num2str(pulseWidth), 'ms)']}];
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

x = [nL4, nL5, nL6, nOther];
explode = ones(1,4);
labels = [{'Put L4'}, {'L5'}, {'Put L6'}, {'Other'}];

figure('Color', 'White');
pie(x, explode)
leg = legend;
leg.String = labels;
leg.Box = 'Off';
leg.Location = 'eastoutside';
leg.FontSize = 8;


ax = gca;
ax.Children(2).FaceColor = [1,1,1];
ax.Children(4).FaceColor = blue;
ax.Children(6).FaceColor = red;
ax.Children(8).FaceColor = green;
ax.FontName = 'Arial';

for slice = 1:2:2*length(x)
ax.Children(slice).FontName = 'Arial';
ax.Children(1).FontSize = 10;
end

ax.Children(2).FaceColor = [1,1,1];
ax.Children(4).FaceColor = blue;
ax.Children(6).FaceColor = red;
ax.Children(8).FaceColor = green;

for slice = 2:2:2*length(x)
    ax.Children(slice).EdgeColor = [1,1,1];
end
%% GroupIDs

GroupIDs = [{'L4'}, ID(L4); {'L5'}, ID(L5); {'L6'}, ID(L6)];
% nGrps = length(GroupIDs);
% for grp = 1:nGrps
   