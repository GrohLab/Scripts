%% Prep Variables
mrControl = clInfo.Mech_Control_4mW_Mech_Rate_Evoked-clInfo.Mech_Control_4mW_preMech_Rate_Evoked;
mrLaser = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked-clInfo.Mech_Laser_5sec_4mW_preMech_Rate_Evoked;
data = table(clInfo.id, clInfo.abs_depth, mrControl, mrLaser);
mrInd = ismember(data.Var1, clInfo.id(clInfo.Mech_Control_4mW_R==true));
data = data(find(mrInd),:);
data = sortrows(data,'Var2','ascend');
ID = data.Var1;
depths = data.Var2;
c = data.Var3;
l = data.Var4;
jitter = randi(20,size(depths));
depths = -1*(depths + jitter);


%% Depth vs Evoked Rate (Control)


deltaRate = zeros(length(depths), 2);
deltaRate(:,2) = c;

figure('Color', 'White', 'Name', 'Changes in FR with Mechanical Stimulation');
hold on
for unit = 1:length(ID)
dr = deltaRate(unit,:);
dpth = [depths(unit), depths(unit)];
if ismember(ID{unit},L6)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 0.5, 1, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Later)
plot(dr, dpth, 'LineStyle', '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Latest)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', 10);
end
end
ax = gca;
ax.Title.String = 'Changes in FR with Mechanical Stimulation';
ax.YLabel.String = 'Depth [\mum]';
ax.XLabel.String = '\DeltaRate [Hz]';
ax.FontName = 'Arial';
ax.FontSize = 30;
ax.YLim = [-1600, 0];
% ax.XLim = [-50, 50];
% ax.XTick = [-50:10:50];
plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');





%% Depth vs Evoked Rate (Laser)


deltaRate = zeros(length(depths), 2);
deltaRate(:,2) = l;

figure('Color', 'White', 'Name', 'Changes in FR with Mechanical + L6 Stimulation');
hold on
for unit = 1:length(ID)
dr = deltaRate(unit,:);
dpth = [depths(unit), depths(unit)];
if ismember(ID{unit},L6)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 0.5, 1, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Later)
plot(dr, dpth, 'LineStyle', '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Latest)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', 10);
end
end
ax = gca;
ax.Title.String = 'Changes in FR with Mechanical + L6 Stimulation';
ax.YLabel.String = 'Depth [\mum]';
ax.XLabel.String = '\DeltaRate [Hz]';
ax.FontName = 'Arial';
ax.FontSize = 30;
ax.YLim = [-1600, 0];
% ax.XLim = [-50, 50];
% ax.XTick = [-50:10:50];
plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');








%%
% Depth vs Delta Rate 8.10.21

deltaRate = zeros(length(depths), 2);
deltaRate(:,2) = l - c;




figure('Color', 'White', 'Name', 'Laser Modulation of Mechanical Responses');

hold on
for unit = 1:length(ID)
dr = deltaRate(unit,:);
dpth = [depths(unit), depths(unit)];
if ismember(ID{unit},L6)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 0.5, 1, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Later)
plot(dr, dpth, 'LineStyle', '-', 'Color', [1, 0, 0, 0.5], 'LineWidth', 10);
elseif ismember(ID{unit}, Latest)
plot(dr, dpth, 'LineStyle', '-', 'Color', [0, 1, 0, 0.5], 'LineWidth', 10);
end
end
ax = gca;
ax.Title.String = 'Laser Modulation of Mechanical Responses';
ax.YLabel.String = 'Depth [\mum]';
ax.XLabel.String = '\DeltaRate [Hz]';
ax.FontName = 'Arial';
ax.FontSize = 30;
ax.YLim = [-1600, 0];
% ax.XLim = [-50, 50];
% ax.XTick = [-50:10:50];
plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');
 