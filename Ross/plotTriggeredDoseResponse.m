function plotTriggeredDoseResponse(ID, sortedData, condIndices, Conditions, fs)

if class(ID) == 'cell'
    figureName = ID{1};
    lngth = 1; % generlaise this!
else
    figureName = ID;
    lngth = 1;
end
rates = zeros(lngth, length(condIndices));
sds = zeros(lngth, length(condIndices));
timeBin = 0.05;

for b= 1:lngth
    spkInd = find(ismember(sortedData(:,1), ID));
    train = sortedData(spkInd, 2);
    for a = 1:length(condIndices)
        Cond{a} = Conditions(condIndices(a)).Triggers;
        rate = TriggeredRates(train, Cond{a},fs, timeBin);
        rates(b,a) = mean(rate{:});
        sds(b,a) = std(rate{:});
    end
end

ax = gca;
errorbar(ax, rates, sds, 'vertical')
xlabel('Conditions');
ax.XLim(2) = length(condIndices) + 1;
ax.XTick = [1:1:length(condIndices)];
xnames = cell(1, length(condIndices));
for a = 1:length(condIndices)
    undInds = strfind(Conditions(condIndices(a)).name, '_');
    Conditions(condIndices(a)).name(undInds) = ' ';
    Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(undInds(1)+1:undInds(end)-1);
    if strfind(Conditions(condIndices(a)).name, '5sec')
        Conditions(condIndices(a)).name = Conditions(condIndices(a)).name(6:end);
        Conditions(condIndices(a)).name(1) = 'P';
    end
    xnames{a} = Conditions(condIndices(a)).name;
end

ax.XTickLabel = xnames;
ax.XTickLabelRotation = -45;
ylabel('Firing Rate _{(Hz)}');
ax.FontSize = 20;
ax.FontName = 'Arial';
end