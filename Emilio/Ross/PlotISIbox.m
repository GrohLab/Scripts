%% Considered conditions selection
% Choose the conditions to look at
% auxSubs = setdiff(1:numel(CFA_ISI.CondName));
% ccondNames = CFA_ISI.CondName(auxSubs);
% [cchCond, iOk] = listdlg('ListString',ccondNames,'SelectionMode','multiple',...
%     'PromptString',...
%     'Choose the condition(s) to look at:',...
%     'ListSize', [350, numel(condNames)*16]);
% if ~iOk
%     fprintf(1,'Cancelling...\n')
%     return
% end

figure('Color',[1,1,1]);
for figNo = 1:3
    subplot(3,1,figNo)
    plot(CFA_ISI(figNo).Vals(2).bns{1}, CFA_ISI(figNo).Vals(2).cts{1}./sum(CFA_ISI(figNo).Vals(2).cts{1}),'LineWidth',1);
    % hold on
    for a = 2: length(CFA_ISI(figNo).Vals(2).cts)
        plot(CFA_ISI(figNo).Vals(2).bns{a}, CFA_ISI(figNo).Vals(2).cts{a}./sum(CFA_ISI(figNo).Vals(2).cts{a}),'LineWidth',1);
    end
    fig = gcf;
%     ax = fig.Children;
%     % ax.XTickLabel = round(exp(cellfun(@str2double,ax.XTickLabel)) * 1e3);
%     %     xlabel(ax,'Time [ms]'); ylabel(ax,'ISI Probability');
%     grid(ax,'on')
    title(CFA_ISI(figNo).CondName)
    
end


