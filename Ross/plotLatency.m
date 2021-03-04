function fig = plotLatency(UnitID, sortedData, TriggerTimes, samplingFrequency)

% Finds TriggerLatencies for a given unit and plots them in box plot
ind = ismember(sortedData(:,1), UnitID);
Latencies = TriggerLatencies(sortedData(ind,2), TriggerTimes, samplingFrequency);
if class(UnitID) == 'cell'
    figureName = UnitID{1};
else
    figureName = num2str(UnitID);
end
fidelity = (length(Latencies{1})/length(TriggerTimes))*100; % how often it responds at least once in a trigger window
fig = figure('Name', ['Unit ', figureName, ' Triggered Latency', ], 'Color', 'white');
boxplot(Latencies{1});
title(fig.Name, 'Interpreter', 'none');
ylabel('Trigger Latency_{(ms)}')
xlim([0.9, 1.1]);
ylim([0,50]);
yticks([0:5:50]);
x = xlim;
y = ylim;
xlow = x(1);
ylow = y(1);
xhigh = x(2);
yhigh = y(2);
xcoord = xlow + (xhigh - xlow)*0.7;
ycoord = ylow + (yhigh - ylow)*0.9;
text(xcoord, ycoord, [num2str(round(fidelity)), '% fidelity'],'FontSize', 15);
ax = gca;
xticks([]);
ax = gca;
ax.FontSize = 15;

end