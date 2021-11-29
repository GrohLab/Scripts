%FR_Modulation

%% FR Modulation Plots
clInd = find(clInfo.Mech_Control_R);

id_depth_mod = table(clInfo.id(clInd), clInfo.abs_depth(clInd), clInfo.Mech_Control_Rate(clInd),...
    clInfo.Mech_Control_Rate_Baseline(clInd), clInfo.Mech_Laser_Rate(clInd), clInfo.Mech_Laser_Rate_Baseline(clInd));
id_depth_mod = sortrows(id_depth_mod,'Var2','ascend');

ID = id_depth_mod{:,1};

mc = id_depth_mod{:,3};
mcb = id_depth_mod{:,4};
ml = id_depth_mod{:,5};
mlb = id_depth_mod{:,6};

deltaRate = [mc-mcb, ml-mlb];



depth = id_depth_mod.Var2;
jitter = randi(20,size(depth));
depth = depth + jitter; % adding small jitter for visualisation
depth = -1*depth;



for chc = 1:2
    
    fig = figure('Name', 'Depth vs Modulation', 'Color', 'White');
    ax = gca;
    hold on
    for cl = 1:length(ID)
        rt = deltaRate(cl,chc);
        dpth = depth(cl);
        plot([0, rt], [dpth, dpth], 'LineStyle', '-', 'Marker', 'none', 'Color', [0, 0, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5]);
        
    end
    
    %     f = fit(depth, modinc, 'smoothingspline');
    %     plot(f,depth, modinc)
    
end

%     legend on
%plot(deltaRate(:,2), depth, 'LineStyle', '-', 'Marker', 'none', 'Color', [0.5, 1, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01,  'MarkerFaceColor',[0.5,0.5,0.5])
title('Depth vs Modulation of Mechanical Responses', 'Interpreter', 'none')
ylabel( 'Depth [\mum]', 'Interpreter', 'tex');
xlabel('Modulation Index', 'Interpreter', 'tex');
%ax.XLim = [-1.1*max(deltaRate(:,2)), 1.1*max(deltaRate(:,2))];
ax.XLim = [-1, 1];

hold on
plot(zeros(length(mod)), linspace(ax.YLim(1),ax.YLim(2), length(mod)), 'LineStyle', '-', 'Color', 'k');
ax.FontSize = 25;
ax.FontName = 'Arial';
%     ax.Legend.String = {'Late Responders'};
%     savefig(fig, fullfile('Z:\Ross\Experiments\smrxFiles\16.12.20\12.6.21\Figures\', [consCondNames{chc}, '_LateResponders_Depth_vs_Modulation.fig']));


