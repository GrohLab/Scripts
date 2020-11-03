%% Saline vs CFA Spontaneous
for model = 1:2
    if model == 1
        Hist = SalISIhist;
    else
        Hist = CfaISIhist;
    end
    
    IsiR = length(Hist(1).Vals(1).cumISI);
    for a = 2: length(Hist)
        IsiR = [IsiR, length(Hist(a).Vals(1).cumISI)];
    end
    r = sum(IsiR);
    Isi = zeros(r, length(Hist(1).Vals(1).cumISI{1}));
    c = 1;
    for a = 1:length(Hist)
        for b = 1:length(Hist(a).Vals(1).cumISI)
            Isi(c,:) = Hist(a).Vals(1).cumISI{b};
            c = c +1;
        end
    end
    % Getting rid of NaNs
    Isi(length(Isi(:,1)) + 1,:) = NaN;
    IndNaN = length(Isi(:,1));
    for nInd = 1:length(Isi(:,1))
        if isnan(Isi(nInd,:)) == true
            IndNaN = [IndNaN; nInd];
        end
    end
    IndNaN = sort(IndNaN, 'descend');
    Isi(IndNaN,:) =[];
    if model == 1
        SpontSalISI = Isi;
        cumSalspont = sum(Isi)/length(Isi(:,1));
    else
        SpontCfaISI = Isi;
        cumCFAspont = sum(Isi)/length(Isi(:,1));
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
%%
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
%% Stats Test and p-Values on Figure
[txt, star] = findKSsignificance(cumSalspont, cumCFAspont);

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
    IsiR = length(Hist(evInd(1)).Vals(1).cumISI);
    for a = 2:length(evInd)
        IsiR = [IsiR, length(Hist(evInd(a)).Vals(2).cumISI)];
    end
    r = sum(IsiR);
    Isi = zeros(r, length(Hist(evInd(1)).Vals(2).cumISI{1}));
    c = 1;
    for a = 1: length(evInd)
        for b = 1:length(Hist(evInd(a)).Vals(2).cumISI)
            Isi(c,:) = Hist(evInd(a)).Vals(2).cumISI{b};
            c = c + 1;
        end
    end
    % Getting rid of NaNs
    Isi(length(Isi(:,1)) + 1,:) = NaN;
    IndNaN = length(Isi(:,1));
    for nInd = 1:length(Isi(:,1))
        if isnan(Isi(nInd,:)) == true
            IndNaN = [IndNaN; nInd];
        end
    end
    IndNaN = sort(IndNaN, 'descend');
    Isi(IndNaN,:) =[];
    if model == 1
        EvokedSalISI = Isi;
        cumSalevoked = sum(Isi)/length(Isi(:,1));
    else
        EvokedCfaISI = Isi;
        cumCFAevoked = sum(Isi)/length(Isi(:,1));
    end
end

figure('Color',[1,1,1], 'Name', 'Saline & CFA Spontaneous vs Evoked Cumulative Fraction');
subplot(2,1,1)
for model = 1:2
    subplot(2,1,model)
    if model == 1
        Hist = SalISIhist;
        spont = cumSalspont;
        evoked = cumSalevoked;
    else
        Hist = CfaISIhist;
        spont = cumCFAspont;
        evoked = cumCFAevoked;
    end
    % Only if the bin widths are constant!!!
    plot(Hist(1).Vals(1).bns{1}, spont);
    hold on
    plot(Hist(1).Vals(2).bns{1}, evoked);
    ylabel('Cumulative Fraction');
    xlabel('ISI (msecs)');
    %     xlim([log(0.001), 4.5]);
    %     ylim([0, 1]);
    legend('Spontaneous', 'Evoked');
    fig = gcf;
    ax = gca;
    ax.FontSize = 20;
    if model == 1
        title('Saline');
    else
        title('CFA');
    end
    %%
    ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
    %% More Stats
    [txt, star] = findKSsignificance(spont, evoked);
    annotation(figure(gcf),'textbox',...
        [0.21428125, 0.438809261300992*model 0.0568125 0.029768467475193],'String',star,...
        'FontSize',20,...
        'FitBoxToText','off',...
        'EdgeColor','none');
    annotation(figure(gcf),'textbox',...
        [0.673958333333333, 0.112328008579419*model 0.203645833333333 0.232874455900804],...
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
            subplot(3,1,a)
            hold on
            for chFig = chCond
                Isi = zeros(length(Histy(chFig).Vals(wIndex).cumISI), length(Histy(chFig).Vals(wIndex).cumISI{1}));
                for cInd = 1:length(Histy(chFig).Vals(wIndex).cumISI)
                    Isi(cInd,:) =  Histy(chFig).Vals(wIndex).cumISI{cInd};
                end
                % Getting rid of NaNs
                Isi(length(Histy(chFig).Vals(wIndex).cumISI) + 1,:) = NaN;
                IndNaN = (length(Histy(chFig).Vals(wIndex).cumISI) + 1);
                for nInd = 1:(length(Histy(chFig).Vals(wIndex).cumISI) + 1)
                    if isnan(Isi(nInd,:)) == true
                        IndNaN = [IndNaN; nInd];
                    end
                end
                IndNaN = sort(IndNaN, 'descend');
                Isi(IndNaN,:) =[];
                cumISI = sum(Isi)/length(Histy(chFig).Vals(wIndex).cumISI);
                % Only if the bin widths are constant!!!
                plot(CfaISIhist(1).Vals(1).bns{1}, cumISI)
                if model == 1
                    SalISIhist(chFig).Vals(wIndex).cumsum = cumISI;
                    Histy = SalISIhist;
                elseif model == 2
                    CfaISIhist(chFig).Vals(wIndex).cumsum = cumISI;
                    Histy = CfaISIhist;
                end
            end
            % Results(c).name = kstest2(cumISI(chCond(1)).Condition(wIndex).Vals;
            ylabel('Cumulative Fraction');
            xlabel('ISI (msecs)');
            xlim([log(0.001), 4.5]);
            ylim([0, 1]);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            fig = gcf;
            ax = gca;
            ax.FontSize = 20;
            ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
           
            
            for a = 1:length(chCond)-1
                for b  = a+1: length(chCond)
                    Results(r).name = [modName, ' ', Histy(a + c).name, ' vs ', Histy(b + c).name, ' ', wName];
                    Results(r).kstest2 = findKSsignificance(Histy(a + c).Vals(wIndex).cumsum,Histy(b + c).Vals(wIndex).cumsum);
                    r = r + 1;
                end
            end
            
            
          c = c + 3;  
        end
        
    end
end
for a = 1:length(Results)
fprintf(1, [Results(a).name,' Kolgorov-Smirnov Test ', num2str(Results(a).kstest2), ' \n \n'])
end