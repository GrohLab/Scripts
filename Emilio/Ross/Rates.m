% Rates
% Spont vs Evoked
SpontInd = find(contains(Vpl.Properties.VariableNames, 'Counts_Spont'));
MechControlInd = contains(Vpl.Properties.VariableNames, 'Mech_Control');
CountsInd = contains(Vpl.Properties.VariableNames, 'Counts_Evoked');
EvInd = find(CountsInd & MechControlInd);
figure('Name', 'CFA Spontaneous vs Evoked Rates (MR Units)', 'Color', 'white');
boxplot([CFA_Spont,  CFA_Evoked], 'Notch', 'on', 'Colors', 'r');
hold on
for i = 1:length(CFA_Spont)
    plot([1.25, 1.75],[CFA_Spont(i,1), CFA_Evoked(i,1)],'-o', 'color', [0.9,0.9,0.9]);
end
ylabel('Firing Rate (Hz)');
labels = [{'Spontaneous'}, {'Evoked'}];
xticklabels(labels);
title('CFA');
[txt, star] = findRSsignificance(CFA_Spont, CFA_Evoked);
annotation(figure(gcf),'textbox',...
    [0.21428125 0.438809261300992 0.0568125 0.029768467475193],'String',star,...
    'FontSize',20,...
    'FitBoxToText','off',...
    'EdgeColor','none');
annotation(figure(gcf),'textbox',...
    [0.673958333333333 0.112328008579419 0.203645833333333 0.232874455900804],...
    'String',{txt},...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');
ax = gca;
ax.FontSize = 25;