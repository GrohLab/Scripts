%% 




%% ISIs
ind = 1;
for model = 1:2
    if model == 1
        Histy = SalIIShist;
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
        
        figure('Color',[1,1,1], 'Name', [modName, ' ', wName, ' ISIs']);
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
            plot(Histy(chCond(1)).Vals(wIndex).bns, Histy(chCond(1)).Vals(wIndex).cts./sum(Histy(chCond(1)).Vals(wIndex).cts),'LineWidth',1);
            hold on
            plot(Histy(chCond(2)).Vals(wIndex).bns, Histy(chCond(2)).Vals(wIndex).cts./sum(Histy(chCond(2)).Vals(wIndex).cts),'LineWidth',1);
            plot(Histy(chCond(3)).Vals(wIndex).bns, Histy(chCond(3)).Vals(wIndex).cts./sum(Histy(chCond(3)).Vals(wIndex).cts),'LineWidth',1);
            ylabel('ISI Probability');
            xlabel('ISI (msecs)');
            % title([num2str(pwr), ' mW']);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            ylim([0, 10^(-2)]);
            xlim([log(0.001), 4.5]);
            fig = gcf;
            ax = gca;
            ax.FontSize = 20;
            ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3); % This goes weird when image is not full size
            
            ISIvals(ind).name = [modName, ' ', wName, ' ', Histy(chCond(1)).name];
            ISIvals(ind).Vals = Histy(chCond(1)).Vals(wIndex).cts./sum(Histy(chCond(1)).Vals(wIndex).cts);
            ind = ind + 1;
            ISIvals(ind).name = [modName, ' ', wName, ' ', Histy(chCond(2)).name];
            ISIvals(ind).Vals = Histy(chCond(2)).Vals(wIndex).cts./sum(Histy(chCond(2)).Vals(wIndex).cts);
            ind = ind + 1;
            ISIvals(ind).name = [modName, ' ', wName, ' ', Histy(chCond(3)).name];
            ISIvals(ind).Vals = Histy(chCond(3)).Vals(wIndex).cts./sum(Histy(chCond(3)).Vals(wIndex).cts);
            ind = ind + 1;
        end
    end
end
%% kstest2
ind = 1;
BonfAlpha = 0.05/length(ISIvals);
c = 0;
for comparisons = 1:3:length(ISIvals)
    for a = 1:2
        for b = a+1:3
    Results(ind).name = [ISIvals(c + a).name, ' vs ', ISIvals(c + b).name,];
    Results(ind).kstest = kstest2(ISIvals(c + a).Vals, ISIvals(c + b).Vals,'Alpha',BonfAlpha);
    ind = ind + 1;
        end
    end
    c = c + 3;
end
for a = 1:36
fprintf(1, [Results(a).name,' Kolgorov-Smirnov Test = ', num2str(Results(a).kstest), ' \n \n'])
end
save(fullfile('Z:\Ross\Experiments\VPL', 'Laser_Results'), 'Results', '-v7.3');
%% Cumulative Fraction
for model = 1:2
    if model == 1
        Histy = SalIsi;
        modName = 'Saline';
    else
        Histy = CFAisi;
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
            Ncts = Histy(chCond(1)).Vals(wIndex).cts/sum(Histy(chCond(1)).Vals(wIndex).cts);
            plot(Histy(chCond(1)).Vals(wIndex).bns,cumsum(Ncts),'LineStyle','-')
            hold on
            Ncts = Histy(chCond(2)).Vals(wIndex).cts/sum(Histy(chCond(2)).Vals(wIndex).cts);
            plot(Histy(chCond(2)).Vals(wIndex).bns,cumsum(Ncts),'LineStyle','-')
            Ncts = Histy(chCond(3)).Vals(wIndex).cts/sum(Histy(chCond(3)).Vals(wIndex).cts);
            plot(Histy(chCond(3)).Vals(wIndex).bns,cumsum(Ncts),'LineStyle','-')
            ylabel('Cumulative Fraction');
            xlabel('ISI (msecs)');
            xlim([log(0.001), 4.5]);
            ylim([0, 1]);
            
            % title([num2str(pwr), ' mW']);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            %             ylim([0, 1]);
            %             xlim([log(0.001), 4.5]);
            ax = gca;
            ax.FontSize = 20;
            grid(ax,'on')
            ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
        end
    end
end
%% Saline vs CFA of Mechanically Responsive Units

Sal.SpontCts = SalIsi(1).Vals(1).cts;
for a = 2:length(SalIsi)
    Sal.SpontCts = Sal.SpontCts + SalIsi(a).Vals(1).cts;
end
Sal.EvokedCts = SalIsi(1).Vals(2).cts;
for a = [4, 7]
    Sal.EvokedCts = Sal.EvokedCts + SalIsi(a).Vals(2).cts;
end
Sal.bns = SalIsi(1).Vals(2).bns;

CFA.SpontCts = CFAisi(1).Vals(1).cts;
for a = 2:length(CFAisi)
    CFA.SpontCts = CFA.SpontCts + CFAisi(a).Vals(1).cts;
end
CFA.EvokedCts = CFAisi(1).Vals(2).cts;
for a = [4, 7]
    CFA.EvokedCts = CFA.EvokedCts + CFAisi(a).Vals(2).cts;
end
CFA.bns = CFAisi(1).Vals(2).bns;

figure('Color',[1,1,1], 'Name', 'Saline vs CFA ISIs');
subplot(2,1,1)
plot(CFA.bns, CFA.SpontCts./sum(CFA.SpontCts),'LineWidth',1);
hold on
plot(Sal.bns, Sal.SpontCts./sum(Sal.SpontCts),'LineWidth',1);

title('Spontaneous');
ylabel('ISI Probability');
xlabel('ISI (msecs)');
legend('CFA', 'Saline');
ylim([0, 10^(-2)]);
xlim([log(0.001), 4.5]);
ax = gca;
ax.FontSize = 20;
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);

subplot(2,1,2)
plot(CFA.bns, CFA.EvokedCts./sum(CFA.EvokedCts),'LineWidth',1);
hold on
plot(Sal.bns, Sal.EvokedCts./sum(Sal.EvokedCts),'LineWidth',1);

title('Evoked');
ylabel('ISI Probability');
xlabel('ISI (msecs)');
legend('CFA', 'Saline');
ylim([0, 10^(-2)]);
xlim([log(0.001), 4.5]);
ax = gca;
ax.FontSize = 20;
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
figure('Color',[1,1,1], 'Name', 'Saline vs CFA Cumulative Fraction');
subplot(2,1,1)
Ncts = CFA.SpontCts/sum(CFA.SpontCts);
plot(CFA.bns,cumsum(Ncts),'LineStyle','-')
hold on
Ncts = Sal.SpontCts/sum(Sal.SpontCts);
plot(Sal.bns,cumsum(Ncts),'LineStyle','-')
title('Spontaneous');
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([log(0.001), 4.5]);
ylim([0, 1]);
legend('CFA', 'Saline');
ax = gca;
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);ax.FontSize = 20;
grid(ax,'on')
subplot(2,1,2)
Ncts = CFA.EvokedCts/sum(CFA.EvokedCts);
plot(CFA.bns,cumsum(Ncts),'LineStyle','-')
hold on
Ncts = Sal.EvokedCts/sum(Sal.EvokedCts);
plot(Sal.bns,cumsum(Ncts),'LineStyle','-')
title('Evoked');
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([log(0.001), 4.5]);
ylim([0, 1]);
legend('CFA', 'Saline');
ax = gca;
ax.FontSize = 20;
ax.XTickLabels = round(exp(cellfun(@str2double,ax.XTickLabel))* 1e3);
grid(ax,'on')
%% Results 
if kstest2(CFA.SpontCts./sum(CFA.SpontCts),Sal.SpontCts./sum(Sal.SpontCts))
    fprintf(1,'Spontaneous activity between Saline and CFA is significantly different (kstest2) \n');
else
    fprintf(1,'Spontaneous activity between Saline and CFA is NOT significantly different (kstest2) \n');
end
%% Laser Differences in Saline vs CFA (Needs to be normalise somehow)



