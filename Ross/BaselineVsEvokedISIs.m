function fig = BaselineVsEvokedISIs(SalISIhist, CfaISIhist)
%% Saline & CFA Spont vs Evoked (Mech Control)
for model = 1:2
    if model == 1
        Hist = SalISIhist;
    else
        Hist = CfaISIhist;
    end
    evInd =[];
    for e = 1:length(Hist)
        if contains(Hist(e).name, 'Mech_Control', 'IgnoreCase', true) == true
            evInd = [evInd; e];
        end
    end
    IsiR = length(Hist(evInd(1)).Vals(1).CumISI);
    for a = 2:length(evInd)
        IsiR = [IsiR, length(Hist(evInd(a)).Vals(2).CumISI)];
    end
    r = sum(IsiR);
    SpontIsiStack = zeros(r, length(Hist(evInd(1)).Vals(1).ISI{1}));
    SpontcumIsiStack = zeros(r, length(Hist(evInd(1)).Vals(1).CumISI{1}));
    EvokedIsiStack = zeros(r, length(Hist(evInd(1)).Vals(2).ISI{1}));
    EvokedcumIsiStack = zeros(r, length(Hist(evInd(1)).Vals(2).CumISI{1}));
    c = 1;
    for a = 1: length(evInd)
        for b = 1:length(Hist(evInd(a)).Vals(2).CumISI)
            SpontcumIsiStack(c,:) = Hist(evInd(a)).Vals(1).CumISI{b};
            EvokedcumIsiStack(c,:) = Hist(evInd(a)).Vals(2).CumISI{b};
            SpontIsiStack(c,:) = Hist(evInd(a)).Vals(1).ISI{b};
            EvokedIsiStack(c,:) = Hist(evInd(a)).Vals(2).ISI{b};
            c = c + 1;
        end
    end
    % Getting rid of NaNs
    SpontcumIsiStack(length(SpontcumIsiStack(:,1)) + 1,:) = NaN;
    SpontIsiStack(length(SpontIsiStack(:,1)) + 1,:) = NaN;
    EvokedcumIsiStack(length(EvokedcumIsiStack(:,1)) + 1,:) = NaN;
    EvokedIsiStack(length(EvokedIsiStack(:,1)) + 1,:) = NaN;
    
    SpontIndNaN = length(SpontIsiStack(:,1));
    EvokedIndNaN = length(EvokedIsiStack(:,1));
    
    
    for nInd = 1:length(SpontIsiStack(:,1))-1
        if isnan(SpontIsiStack(nInd,:)) == true
            SpontIndNaN = [SpontIndNaN; nInd];
        end
    end
    
    for nInd = 1:length(EvokedIsiStack(:,1))-1
        if isnan(EvokedIsiStack(nInd,:)) == true
            EvokedIndNaN = [EvokedIndNaN; nInd];
        end
    end
    SpontIndNaN = sort(SpontIndNaN, 'descend');
    SpontIsiStack(SpontIndNaN,:) =[];
    SpontcumIsiStack(SpontIndNaN,:) = [];
    EvokedIndNaN = sort(EvokedIndNaN, 'descend');
    EvokedIsiStack(EvokedIndNaN,:) =[];
    EvokedcumIsiStack(EvokedIndNaN,:) = [];
    if model == 1
        cumSpontSalIsi = sum(SpontcumIsiStack)/length(SpontcumIsiStack(:,1));
        cumEvokedSalIsi = sum(EvokedcumIsiStack)/length(EvokedcumIsiStack(:,1));
        sumSpontSalISI = sum(SpontIsiStack)/length(SpontIsiStack(:,1));
        sumEvokedSalISI = sum(EvokedIsiStack)/length(EvokedIsiStack(:,1));
    else
        cumSpontCfaIsi = sum(SpontcumIsiStack)/length(SpontcumIsiStack(:,1));
        cumEvokedCfaIsi = sum(EvokedcumIsiStack)/length(EvokedcumIsiStack(:,1));
        sumSpontCfaISI = sum(SpontIsiStack)/length(SpontIsiStack(:,1));
        sumEvokedCfaISI = sum(EvokedIsiStack)/length(EvokedIsiStack(:,1));
    end
end

figure('Color',[1,1,1], 'Name', 'Saline & CFA Spontaneous vs Evoked Cumulative Fraction');
subplot(2,1,1)
for model = 1:2
    subplot(2,1,model)
    if model == 1
        Hist = SalISIhist;
        spont = cumSpontSalIsi;
        evoked = cumEvokedSalIsi;
        SpontISI = sumSpontSalISI;
        EvokedISI = sumEvokedSalISI;
    else
        Hist = CfaISIhist;
        spont = cumSpontCfaIsi;
        evoked = cumEvokedCfaIsi;
        SpontISI = sumSpontCfaISI;
        EvokedISI = sumEvokedCfaISI;
    end
    % Only if the bin widths are constant!!!
    plot(Hist(1).Vals(1).bns{1}, spont);
    hold on
    plot(Hist(1).Vals(2).bns{1}, evoked);
    ylabel('Cumulative Fraction');
    xlabel('ISI (msecs)');
    xlim([log(0.001), 4.5]);
    ylim([0, 1]);
    legend('Spontaneous', 'Evoked');
    fig = gcf;
    ax = gca;
    ax.FontSize = 20;
    if model == 1
        title('Saline');
    else
        title('CFA');
    end
    ax.XTickLabels = round(exp(-6:4)* 1e3);
    %% More Stats
    [txt, star] = findKSsignificance(SpontISI, EvokedISI);
    annotation(figure(gcf),'textbox',...
        [0.21428125, 0.438809261300992*(1/model) 0.0568125 0.029768467475193],'String',star,...
        'FontSize',20,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    annotation(figure(gcf),'textbox',...
        [0.673958333333333, 0.112328008579419 *(1/model) 0.203645833333333 0.232874455900804],...
        'String',{txt},...
        'FontSize',15,...
        'FitBoxToText','off',...
        'EdgeColor','none');
end
end