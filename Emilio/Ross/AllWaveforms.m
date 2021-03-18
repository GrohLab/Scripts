function fig = AllWaveforms(clWaveforms, ID, fs)
fig = figure('color', 'white');
hold on
for a = 1:length(ID)
    ind = find(ismember(clWaveforms(:,1), ID(a)));
    plot(mean(clWaveforms{ind,2}'));
end
fig = gcf;
ax = fig(1).Children;
ind = size(clWaveforms{1,2}(:,1));
ind = ind(1);
ax.XTickLabel = [1:ind]/fs*1000;
title('m379 All Units')
xlabel(ax,'Time_(_m_s_)')
ax.FontSize = 30;
end