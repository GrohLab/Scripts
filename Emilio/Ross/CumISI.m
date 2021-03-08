

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
figure('Color',[1,1,1], 'Name', 'Saline vs CFA Spontaneous Cumulative Fraction');
% Only if the bin widths are constant!!!
plot(SalISIhist(1).Vals(1).bns{1}, cumSalspont);
hold on
plot(CfaISIhist(1).Vals(1).bns{1}, cumCFAspont);
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([log(0.001), 4.5]);
ylim([0, 1]);
legend('Saline', 'CFA');
title('Spontaneous Activity Cumulative Fractions');
fig = gcf;
ax = gca;
ax.FontSize = 20;

ax.XTickLabels = round(exp(-6:4)* 1e3);
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

%% CumISI Comparing Conditions
r = 1;
for model = 1:2
    if model == 1
        Histy = SalISIhist;
        modName = 'Saline';
    else
        Histy = CfaISIhist;
        modName = 'CFA';
    end
    for wIndex = 1:2
        if wIndex == 1
            wName = 'Spontaneous';
        else
            wName = 'Evoked';
        end
        
        figure('Color',[1,1,1], 'Name', [modName, ' ', wName, ' Cumulative Fraction']);
        
        
        for a = 1:3
            pwr = a*5;
            if a == 1
                chCond = [7 8 9];
            elseif a == 2
                chCond = [1 2 3];
            elseif a == 3
                chCond = [4 5 6];
            end
            subplot(3,1,a)
            hold on
            for chFig = chCond
                IsiStack = zeros(length(Histy(chFig).Vals(wIndex).ISI), length(Histy(chFig).Vals(wIndex).ISI{1}));
                cumIsiStack = zeros(length(Histy(chFig).Vals(wIndex).CumISI), length(Histy(chFig).Vals(wIndex).CumISI{1}));
                for cInd = 1:length(Histy(chFig).Vals(wIndex).CumISI)
                    IsiStack(cInd,:) =  Histy(chFig).Vals(wIndex).ISI{cInd};
                    cumIsiStack(cInd,:) =  Histy(chFig).Vals(wIndex).CumISI{cInd};
                end
                % Getting rid of NaNs
                cumIsiStack(length(Histy(chFig).Vals(wIndex).CumISI) + 1,:) = NaN;
                IsiStack(length(Histy(chFig).Vals(wIndex).ISI) + 1,:) = NaN;
                cumIndNaN = (length(Histy(chFig).Vals(wIndex).CumISI) + 1);
                for nInd = 1:(length(Histy(chFig).Vals(wIndex).CumISI))
                    if isnan(cumIsiStack(nInd,:)) == true
                        cumIndNaN = [cumIndNaN; nInd];
                    end
                end
                cumIndNaN = sort(cumIndNaN, 'descend');
                cumIsiStack(cumIndNaN,:) =[];
                IsiStack(cumIndNaN,:) = [];
                ISIcum = sum(cumIsiStack)/length(Histy(chFig).Vals(wIndex).CumISI);
                % Only if the bin widths are constant!!!
                plot(CfaISIhist(1).Vals(1).bns{1}, ISIcum)
                if model == 1
                    SalISIhist(chFig).Vals(wIndex).cumsum = ISIcum;
                    SalISIhist(chFig).Vals(wIndex).ISIsum = IsiStack/length(Histy(chFig).Vals(wIndex).ISI);
                    % SalISIhist(chFig).Vals(wIndex).sumISI = sum(Isi);
                    Histy = SalISIhist;
                elseif model == 2
                    CfaISIhist(chFig).Vals(wIndex).cumsum = ISIcum;
                    CfaISIhist(chFig).Vals(wIndex).ISIsum = IsiStack/length(Histy(chFig).Vals(wIndex).ISI);
                    % CfaISIhist(chFig).Vals(wIndex).sumISI = sum(Isi);
                    Histy = CfaISIhist;
                end
            end
            % Results(c).name = kstest2(CumISI(chCond(1)).Condition(wIndex).Vals;
            ylabel('Cumulative Fraction');
            xlabel('ISI (msecs)');
            xlim([log(0.001), 4.5]);
            ylim([0, 1]);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            fig = gcf;
            ax = gca;
            ax.FontSize = 20;
            ax.XTickLabels = round(exp(-6:4)* 1e3);

        end
        c = 0;
        for a = 1:3
            pwr = a*5;
            if a == 1
                chCond = [7 8 9];
            elseif a == 2
                chCond = [1 2 3];
            elseif a == 3
                chCond = [4 5 6];
            end
            for d = 1:length(chCond)-1
                for b  = d+1: length(chCond)
                    ksResults(r).name = [modName, ' ', Histy(chCond(d) + c).name, ' vs ', Histy(chCond(b) + c).name, ' ', wName];
                    ksResults(r).kstest2 = findKSsignificance(sum(Histy(chCond(d) + c).Vals(wIndex).ISIsum), sum(Histy(chCond(b) + c).Vals(wIndex).ISIsum));
                    r = r + 1;
                end
            end
          
        end
    end
end
%%
for a = 1:length(ksResults)
    fprintf(1, [ksResults(a).name,' Kolgorov-Smirnov Test ', num2str(ksResults(a).kstest2), ' \n \n'])
end