chCond = 1;
figure('Color', 'White', 'Name', [TrigISIs(chCond).name]);
hold on
for wIndex = 1:2
    CondNames = {'Baseline', 'Evoked'};
    clr = {[0.5, 0, 0], [0, 1, 1]};
    IsiStack = zeros(length(TrigISIs(chCond).Vals(wIndex).ISI), length(TrigISIs(chCond).Vals(wIndex).ISI{1}));
    cumIsiStack = zeros(length(TrigISIs(chCond).Vals(wIndex).CumISI), length(TrigISIs(chCond).Vals(wIndex).CumISI{1}));
    for cInd = 1:length(TrigISIs(chCond).Vals(wIndex).CumISI)
        IsiStack(cInd,:) =  TrigISIs(chCond).Vals(wIndex).ISI{cInd};
        cumIsiStack(cInd,:) =  TrigISIs(chCond).Vals(wIndex).CumISI{cInd};
    end
    % Getting rid of NaNs
    cumIsiStack(length(TrigISIs(chCond).Vals(wIndex).CumISI) + 1,:) = NaN;
    IsiStack(length(TrigISIs(chCond).Vals(wIndex).ISI) + 1,:) = NaN;
    cumIndNaN = (length(TrigISIs(chCond).Vals(wIndex).CumISI) + 1);
    for nInd = 1:(length(TrigISIs(chCond).Vals(wIndex).CumISI))
        if isnan(cumIsiStack(nInd,:)) == true
            cumIndNaN = [cumIndNaN; nInd];
        end
    end
    cumIndNaN = sort(cumIndNaN, 'descend');
    cumIsiStack(cumIndNaN,:) =[];
    IsiStack(cumIndNaN,:) = [];
    stckSz = size(cumIsiStack);
    ISIcum = sum(cumIsiStack)/stckSz(1);
    % Only if the bin widths are constant!!!
    plot(TrigISIs(1).Vals(1).bns{1}, ISIcum, 'Color', clr{wIndex})
end
legend(CondNames)
ylabel('Cumulative Fraction');
xlabel('ISI (msecs)');
xlim([-3, 1]);
xticks([-3:1])
ylim([0, 1]);
fig = gcf;
ax = gca;
ax.FontSize = 20;
ax.XTickLabel = 10.^cellfun(@str2double,ax.XTickLabel) * 1e3;