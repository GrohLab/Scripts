function [fig, ksResults] = CondISIs(SalISIhist, CfaISIhist)
fig = gobjects(4);
%% CumISI Comparing Conditions
r = 1;
cf = 1;
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
        
        fig(cf) = figure('Color',[1,1,1], 'Name', [modName, ' ', wName, ' Cumulative Fraction']);
        cf = cf + 1;
        
        for a = 3%:3
            pwr = a*5;
            if a == 1
                chCond = [7 8 9];
            elseif a == 2
                chCond = [1 2 3];
            elseif a == 3
                chCond = [4 5 6];
            end
            % subplot(3,1,a)
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
            xlim([-3, 1]);
            xticks([-3:1])
            ylim([0, 1]);
            legend(Histy(chCond(1)).name, Histy(chCond(2)).name, Histy(chCond(3)).name);
            fig = gcf;
            ax = gca;
            ax.FontSize = 20;
            ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;

        end
        c = 0;
%         for a = 1:3
%             pwr = a*5;
%             if a == 1
%                 chCond = [7 8 9];
%             elseif a == 2
%                 chCond = [1 2 3];
%             elseif a == 3
%                 chCond = [4 5 6];
%             end
%             for d = 1:length(chCond)-1
%                 for b  = d+1: length(chCond)
%                     ksResults(r).name = [modName, ' ', Histy(chCond(d) + c).name, ' vs ', Histy(chCond(b) + c).name, ' ', wName];
%                     ksResults(r).kstest2 = findKSsignificance(sum(Histy(chCond(d) + c).Vals(wIndex).ISIsum), sum(Histy(chCond(b) + c).Vals(wIndex).ISIsum));
%                     r = r + 1;
%                 end
%             end
          
        end
    end
end
%%
% for a = 1:length(ksResults)
%     fprintf(1, [ksResults(a).name,' Kolgorov-Smirnov Test ', num2str(ksResults(a).kstest2), ' \n \n'])
% end