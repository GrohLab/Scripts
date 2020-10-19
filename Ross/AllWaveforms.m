function fig = AllWaveforms(clWaveforms, fs)

fig = figure('color', 'white'); plot(mean(clWaveforms{1,2}'));
hold on
for a = 2:length(clWaveforms)
plot(mean(clWaveforms{a,2}'));
end
fig = gcf;
ax = fig(1).Children;
ind = size(clWaveforms{1,2}(:,1));
ind = ind(1);
ax.XTickLabel = [1:ind]/fs*1000;
title('All Waveforms')
xlabel(ax,'Time [ms]')
end