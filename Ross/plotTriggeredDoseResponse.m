function plotTriggeredDoseResponse(ID, sortedData, condIndices, Conditions, fs)



% atm only works on block conditions....come back to this


if class(ID) == 'cell'
    figureName = ID{1};
    lngth = 1; % generalise this!
else
    figureName = ID;
    lngth = 1;
end
rates = zeros(lngth, length(condIndices));
sds = zeros(lngth, length(condIndices));
timeBin = 5;
conds = zeros(1,length(condIndices));
for b = 1:lngth
    spkInd = ismember(sortedData(:,1), ID);
    train = sortedData(spkInd, 2);
    Cond = cell(1, length(condIndices));
    for a = 1:length(condIndices)
        Cond{a} = Conditions(condIndices(a)).Triggers;
        rate = TriggeredRates(train, Cond{a},fs, timeBin);
        rates(b,a) = mean(rate{:});
        sds(b,a) = std(rate{:});
        if contains(Conditions(condIndices(a)).name, '1Hz') % Can make it so we don't need to know text i.e. Hz Ind = yadayadayada...later
            conds(a) = 1;
        elseif contains(Conditions(condIndices(a)).name, '10Hz')
            conds(a) = 2;
        elseif contains(Conditions(condIndices(a)).name, 'pulse')
            conds(a) = 3;
        end
    end
end

cond1 = conds == 1;
cond2 = conds == 2;
cond3 = conds == 3;
legText = [];
legendCell = [{['1Hz', legText]}, {['10Hz', legText]}, {['Pulse', legText]}];
ax = gca;



xnames = cell(1, length(condIndices));
pwrs = zeros(1, length(condIndices));
for a = 1:length(condIndices)
    undInds = strfind(Conditions(condIndices(a)).name, '_');
    Conditions(condIndices(a)).name(undInds) = ' ';
    Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(undInds(1)+1:undInds(end)-1);
    if strfind(Conditions(condIndices(a)).name, '5sec')
        Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(6:end);
        Conditions(condIndices(a)).name(1) = 'P';
    end
    xnames{a} = Conditions(condIndices(a)).name;
    spaceInd = strfind(xnames{a}, ' ');
    mWind = strfind(xnames{a}, 'mW');
    pwrs(a) = str2double(xnames{a}(spaceInd+1:mWind-1));
end
pwrs = unique(pwrs);
for a = 1:length(pwrs)
pwrCell{a} = num2str(pwrs(a));
end

hold on
x = 1:length(pwrs);
offset = 0.1;
errorbar(ax, x-offset, rates(cond1), sds(cond1), 'vertical', 'LineWidth', 3)
errorbar(ax, x, rates(cond2), sds(cond2), 'vertical', 'LineWidth', 3)
errorbar(ax, x+offset, rates(cond3), sds(cond3), 'vertical', 'LineWidth', 3)

legend(legendCell, 'Location', 'northwest')
xlabel('Laser Intensity _{(mW)}');
ax.XLim(1) = 0;
ax.XLim(2) = length(pwrs) + 1;
ax.XTick = [1:1:length(pwrs)];
ax.XTickLabel = pwrCell;
ylabel('Firing Rate _{(Hz)}');
ax.FontSize = 20;
ax.FontName = 'Arial';
end