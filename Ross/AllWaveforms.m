function fig = AllWaveforms(clWaveforms, ID, fs)
fig = figure('color', 'white'); plot(mean(clWaveforms{ind(1),2}'));
hold on
for a = length(ID)
    ind = find(ismember(clWaveforms(:,1), id(a)));
    plot(mean(clWaveforms{ind,2}'));
end
fig = gcf;
ax = fig(1).Children;
ind = size(clWaveforms{1,2}(:,1));
ind = ind(1);
ax.XTickLabel = [1:ind]/fs*1000;
title('Increased')
xlabel(ax,'Time [ms]')
ax.FontSize = 20;
end