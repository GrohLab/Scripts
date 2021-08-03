%function [increased, decreased, fig] = DepthFrModulation(clInfo, consideredConditions, consCondNames, Counts, IDs, gclID)
clInd = ismember(clInfo.id, nonTagged);
depths = table(clInfo.id(clInd), clInfo.abs_depth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
depth = depths.Var2;
jitter = randi(20,size(depth));
depth = depth + jitter; % adding small jitter for visualisation
depth = -1*depth;
for chc = 15
    spontCounts = Counts{chc,1};
    evCounts = Counts{chc,2};
    c = size(spontCounts);
    c = c(2);
    ordSpontCounts = zeros(length(ID), c);
    ordEvCounts = zeros(length(ID), c);
    for cl = 1:length(ID)
        ind = find(ismember(gclID, ID(cl)));
        ordSpontCounts(cl,:) = spontCounts(ind,:);
        ordEvCounts(cl,:) = evCounts(ind,:);
    end
    spWindow = spontaneousWindow(2)-spontaneousWindow(1);
    evWindow = responseWindow(2)-responseWindow(1);
    medSpontRate = median(ordSpontCounts,2)/spWindow;
    medEvokedRate = median(ordEvCounts,2)/evWindow;
    deltaRate = zeros(length(depth), 2);
    deltaRate(:,2) = medEvokedRate - medSpontRate;
    increased{chc} = ID(deltaRate(:,2)>0);
    decreased{chc} = ID(deltaRate(:,2)<0);
    fig = figure('Name', 'Depth vs Modulation', 'Color', 'White');
    ax = gca;
    hold on
    for cl = 1:length(ID)
        dr = deltaRate(cl,:);
        dpth = [depth(cl), depth(cl)];
        plot(dr, dpth, 'LineStyle', '-', 'Marker', 'none', 'Color', [0, 0.5, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5]);
    end
    
    
    
    %plot(deltaRate(:,2), depth, 'LineStyle', '-', 'Marker', 'none', 'Color', [0.5, 1, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5])
    title(['Depth vs Change in Firing Rate: ', consCondNames{chc}], 'Interpreter', 'none')
    ylabel( 'Depth_{(\mum)}', 'Interpreter', 'tex');
    xlabel('\Delta Firing Rate_{ Hz}', 'Interpreter', 'tex');
    ax.XLim = [-1.1*max(deltaRate(:,2)), 1.1*max(deltaRate(:,2))];
    % ax.YLim = [ax.YLim(1), 0];
    
    hold on
    plot(deltaRate(:,1), linspace(ax.YLim(1),ax.YLim(2), length(deltaRate)), 'LineStyle', '-', 'Color', 'k');
    ax.FontSize = 20;
% end




end
