%% Firing Rates Violin Plots

MechResponse = clInfo.Mech_Control_4mW_R;
ids = [{'L6'}, {'S1'}, {'Vpl'}];
idxMat = [idxTagged, idxNonTagged]; % & MechResponse;
matDims = size(idxMat);
figure('Name', 'Firing Rate Violin Plots', 'Color', 'White');
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s), max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s)];
yMax = max(yMax);
yMin = 0;

for unitGroup = 1:nGroups
    subplot(1, nGroups, unitGroup)
    
    unitIds = idxMat(:,unitGroup);
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s(unitIds);
    rts(:,2) = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s(unitIds);
    
    subplot(1,2,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.2);
       
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
   
    
    plot(zeros(nUnits,1) - 0.01, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', [1,0,0], 'MarkerSize', 5);
    
    plot(zeros(nUnits,1) + 0.01, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 5);
    
    
    
    legend('Mech Control','Mech Laser', 'Location','northwest');
    ax = gca;
    ax.Title.String = ids{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 25;
    ax.YAxis.Limits = [yMin, yMax];
end

%% Relative Rates Violin Plots

MechResponse = clInfo.Mech_Control_4mW_R;
ids = [{'L6'}, {'S1'}, {'Vpl'}];
idxMat = [idxTagged, idxNonTagged]; % & MechResponse;
matDims = size(idxMat);
figure('Name', 'Relative Firing Rate Violin Plots', 'Color', 'White');
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s),...
 max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s)];
yMax = max(yMax);
yMin = [min(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s),...
 min(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s)];
yMin = min(yMin);


for unitGroup = 1:nGroups
    subplot(1, nGroups, unitGroup)
    
    unitIds = idxMat(:,unitGroup);
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Control_4mW_Rate_Spont_1_to_2s(unitIds);
    rts(:,2) = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_1_to_2s(unitIds);
    
    subplot(1,2,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.2);
       
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
   
    
    plot(zeros(nUnits,1) - 0.01, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', [1,0,0], 'MarkerSize', 5);
    
    plot(zeros(nUnits,1) + 0.01, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 5);
    
    
    
    legend('Mech Control','Mech Laser', 'Location','northwest');
    
    ax = gca;
    ax.Title.String = ids{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 25;
    ax.YAxis.Limits = [yMin, yMax];  
    
end
%% SNR Violin Plots


MechResponse = clInfo.Mech_Control_4mW_R;
ids = [{'L6'}, {'S1'}, {'Vpl'}];
idxMat = [idxTagged, idxNonTagged]; % & MechResponse;
matDims = size(idxMat);
figure('Name', 'SNR Violin Plots', 'Color', 'White');
nGroups = matDims(2);

yMax = [max(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s),...
 max(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s)];
yMax = max(yMax);
yMin = [min(clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s),...
 min(clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s)];
yMin = min(yMin);

for unitGroup = 1:nGroups
    subplot(1, nGroups, unitGroup)
    
    unitIds = idxMat(:,unitGroup);
    nUnits = sum(unitIds);
    rts = NaN(nUnits, 2);
    rts(:,1) = clInfo.Mech_Control_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Control_4mW_Rate_Spont_0_to_1s(unitIds);
    rts(:,2) = clInfo.Mech_Laser_5sec_4mW_Rate_Evoked_1_to_2s(unitIds) - clInfo.Mech_Laser_5sec_4mW_Rate_Spont_0_to_1s(unitIds);
    
    subplot(1,2,unitGroup)
    hold on
    
    
    [yCtr, xCtr] = ksdensity(rts(:,1),'Bandwidth',0.7);
    patch(yCtr * -1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.2);
       
    [yCtr, xCtr] = ksdensity(rts(:,2),'Bandwidth',0.7);
    patch(yCtr,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.2);
    
   
    
    plot(zeros(nUnits,1) - 0.01, rts(:,1),...
        'LineStyle','none', 'Marker', '*', 'Color', [1,0,0], 'MarkerSize', 5);
    
    plot(zeros(nUnits,1) + 0.01, rts(:,2),...
        'LineStyle','none', 'Marker', '*', 'Color', [0,0,1], 'MarkerSize', 5);
    
    
    
    legend('Mech Control','Mech Laser', 'Location','northwest');
    
    ax = gca;
    ax.Title.String = ids{unitGroup};
    ax.FontName = 'Arial';
    ax.FontSize = 25;
    ax.YAxis.Limits = [yMin, yMax];  
    
end