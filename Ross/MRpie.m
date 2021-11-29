%% MR Pie Charts

red = [0.75, 0, 0];
green = [0, 0.75, 0];
blue = [0.25, 0.5, 1];

ruNames = {'VPL'}; %[{'Group 1'}, {'Group 2'}, {'L6'}];

% LatertblInd = ismember(clInfo.id, Later);
% LatesttblInd = ismember(clInfo.id, Latest);
% L6tblInd = ismember(clInfo.id, L6);
mrControl = clInfo.Mech_Control_15mW_MR;
mrLaser = clInfo.Mech_Laser_15mW_MR;
mCinc = clInfo.Mech_Control_15mW_Counts_Evoked > clInfo.Mech_Control_15mW_Counts_Spont; mCinc = mCinc/diff(responseWindow);
mLinc = clInfo.Mech_Laser_15mW_Counts_Evoked > clInfo.Mech_Laser_15mW_Counts_Spont; mLinc = mLinc/diff(responseWindow);

inds = clInfo.ActiveUnit; % [LatertblInd, LatesttblInd, L6tblInd];
colours = [red; green; blue];
purple = [0.5,0,0.5];
ruIdxDim = size(ruIdx);
nFilters = ruIdxDim(2);

for p = 1:nFilters
    nCl = sum(inds(:,p));
    nMRcontrol = sum(inds(:,p) & mrControl & mCinc);
    nMRlaser = sum(inds(:,p) & mrLaser & mLinc);
    
    figure('Color', 'White', 'Name', [ruNames{p}, '_MR_PieCharts']);
    
    subplot(2,1,1)
    dataControl = [nCl-nMRcontrol, nMRcontrol];
    explode = [1,1];
    pie(dataControl, explode)
    
    ax = gca;
    axSz = size(ax.Children,1);
    if axSz > 2
        ax.Children(2).EdgeColor = [1,1,1];
        ax.Children(2).FaceColor = purple;
        ax.Children(4).EdgeColor = [1,1,1];
        ax.Children(4).FaceColor = 0.5 * ones(1,3);
        
    else
        ax.Children(2).EdgeColor = [1,1,1];
        ax.Children(2).FaceColor = 0.5 * ones(1,3);
    end
    
    leg = legend({'Non-Responsive', 'Increased Firing Rate [Hz]'});
    leg.FontSize = 10;
    leg.Location = 'eastoutside';
    leg.Box = 'off';
    ax.FontName = 'Arial';
    ax.FontSize = 15;
    ax.Title.String = 'Mech Control';
    
    
    subplot(2,1,2)
    dataLaser = [nCl-nMRlaser, nMRlaser];
    explode = [1,0];
    pie(dataLaser, explode)
    
    ax = gca;
    axSz = size(ax.Children,1);
    if axSz > 2
        ax.Children(2).EdgeColor = [1,1,1];
        ax.Children(2).FaceColor = purple;
        ax.Children(4).EdgeColor = [1,1,1];
        ax.Children(4).FaceColor = 0.5 * ones(1,3);
        
    else
        ax.Children(2).EdgeColor = [1,1,1];
        ax.Children(2).FaceColor = 0.5 * ones(1,3);
    end
    
    
    ax.Title.String = 'Mech + Laser';
    leg = legend({'Non-Responsive', 'Increased Firing Rate [Hz]'});
    leg.Location = 'eastoutside';
    leg.Box = 'off';
    leg.FontSize = 10;
    ax.FontName = 'Arial';
    ax.FontSize = 15;
end