%% Prep Variables
mrControl = clInfo.Mech_Control_4mW_Rate_Evoked-clInfo.Mech_Control_4mW_Rate_Spont;
mrLaser = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked-clInfo.Mech_Laser_5sec_4mW_Rate_Spont;
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

plottedData = [{'Mech'}, {[zeros(length(c),1), c]};...
    {'Mech + L6'}, {[zeros(length(l),1), l]};...
    {'L6 Modulation of Mechanical Responses'}, {[zeros(length(l),1), l-c]}];


%% Plotting

red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];

figure('Color', 'White');
firstLater = find(ismember(ID,Later), 1, 'first');
firstLatest = find(ismember(ID,Latest), 1, 'first');
firstL6 = find(ismember(ID,L6), 1, 'first');
for fg = 1:3
    subplot(2,2,fg);
    if fg == 3
        hold on
        subplot(2,2,3:4)
    end
    hold on
    
    % Plotting the first of each group for the legend
    dr = plottedData{fg,2}(firstLater,:);
    dpth = [depths(firstLater), depths(firstLater)];
    plot(dr, dpth, 'LineStyle', '-', 'Color', [red, 0.5], 'LineWidth', 10);
    
    dr = plottedData{fg,2}(firstLatest,:);
    dpth = [depths(firstLatest), depths(firstLatest)];
    plot(dr, dpth, 'LineStyle', '-', 'Color', [green, 0.5], 'LineWidth', 10);
    
    dr = plottedData{fg,2}(firstL6,:);
    dpth = [depths(firstL6), depths(firstL6)];
    plot(dr, dpth, 'LineStyle', '-', 'Color', [blue, 0.5], 'LineWidth', 10);
    if fg == 3
        leg = legend;
        leg.String = [{'Group 1'}, {'Group 2'}, {'L6'}];
    end
    
    
    
    
    
    
    for unit = 1:length(ID)
        dr = plottedData{fg,2}(unit,:);
        dpth = [depths(unit), depths(unit)];
        
        if ismember(ID{unit}, Later) & unit ~= firstLater
            plot(dr, dpth, 'LineStyle', '-', 'Color', [red, 0.5], 'LineWidth', 10);
        elseif ismember(ID{unit}, Latest) & unit ~= firstLatest
            plot(dr, dpth, 'LineStyle', '-', 'Color',[green, 0.5], 'LineWidth', 10);
        elseif ismember(ID{unit}, L6) & unit ~= firstL6
            plot(dr, dpth, 'LineStyle', '-', 'Color', [blue 0.5], 'LineWidth', 10);
        end
    end
    ax = gca;
    ax.Title.String = plottedData{fg,1};
    ax.YLabel.String = 'Depth [\mum]';
    ax.XLabel.String = '\DeltaRate [Hz]';
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    ax.YLim = [-1600, 0];
    % ax.XLim = [-50, 50];
    % ax.XTick = [-50:10:50];
    plot(zeros(length(ID)), linspace(ax.YLim(1),ax.YLim(2), length(ID)), 'LineStyle', '-', 'Color', 'k');
    if fg == 3
        leg.String(4:end) = [];
        leg.String(4:end) = [];
        leg.Box = 'off';
        legend('Location', 'southeastoutside');
    end
end