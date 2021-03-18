function fig = SalineVsCFAISI(SalISIhist, CfaISIhist)
%% Saline vs CFA Spontaneous
for model = 1:2
    if model == 1
        Hist = SalISIhist;
    else
        Hist = CfaISIhist;
    end
    
    IsiR = length(Hist(1).Vals(1).CumISI);
    for a = 2: length(Hist)
        IsiR = [IsiR, length(Hist(a).Vals(1).CumISI)];
    end
    r = sum(IsiR);
    cumIsiStack = zeros(r, length(Hist(1).Vals(1).CumISI{1}));
    IsiStack = zeros(r, length(Hist(1).Vals(1).ISI{1}));
    c = 1;
    for a = 1:length(Hist)
        for b = 1:length(Hist(a).Vals(1).CumISI)
            cumIsiStack(c,:) = Hist(a).Vals(1).CumISI{b};
            IsiStack(c,:) = Hist(a).Vals(1).ISI{b};
            c = c +1;
        end
    end
    % Getting rid of NaNs
    cumIsiStack(length(cumIsiStack(:,1)) + 1,:) = NaN;
    IsiStack(length(IsiStack(:,1)) + 1,:) = NaN;
    cumIndNaN = length(cumIsiStack(:,1));
    IndNaN = length(IsiStack(:,1));
    for nInd = 1:length(cumIsiStack(:,1))
        if isnan(cumIsiStack(nInd,:)) == true
            cumIndNaN = [cumIndNaN; nInd];
        end
        if isnan(IsiStack(nInd,:)) == true
            IndNaN = [IndNaN; nInd];
        end
    end
    cumIndNaN = sort(cumIndNaN, 'descend');
    cumIsiStack(cumIndNaN,:) =[];
    IndNaN = sort(IndNaN, 'descend');
    IsiStack(IndNaN,:) =[];
    if model == 1
        sumSalIsi = sum(IsiStack)/length(IsiStack(:,1));
        cumSalspont = sum(cumIsiStack)/length(cumIsiStack(:,1));
        
    else
        sumCfaIsi = sum(IsiStack)/length(IsiStack(:,1));
        cumCFAspont = sum(cumIsiStack)/length(cumIsiStack(:,1));
    end
end
fig = figure('Color',[1,1,1], 'Name', 'Saline vs CFA Spontaneous Cumulative Fraction');
% Only if the bin widths are constant!!!
plot(SalISIhist(1).Vals(1).bns{1}, cumSalspont);
hold on
plot(CfaISIhist(1).Vals(1).bns{1}, cumCFAspont);
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([-3, 2]);
xticks([-3:2]);
ylim([0, 1]);
legend('Saline', 'CFA');
title('Spontaneous Activity Cumulative Fractions');
fig = gcf;
ax = gca;
ax.FontSize = 20;
ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;
ax.XTickLabelRotation = 45;
%% Stats Test and p-Values on Figure
[txt, star] = findKSsignificance(sumSalIsi, sumCfaIsi);

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
end
