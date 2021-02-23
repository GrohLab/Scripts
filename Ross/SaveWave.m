function fig = SaveWave(clWaveforms, ID, fs, figureDir)
WaveformDir = fullfile(figureDir,'Waveforms\');
for a = 1:length(ID)
    fig = figure('color', 'white');
    ind = find(ismember(clWaveforms(:,1), ID(a)));
    plot(mean(clWaveforms{ind,2}'));
    fig = gcf;
    ax = fig(1).Children;
    ind = size(clWaveforms{1,2}(:,1));
    ind = ind(1);
    ax.XTick = [0:20:77];
    ax.XTickLabel = [1:ind]/fs*1000;
    title(['Unit ', num2str(ID{a})])
    xlabel(ax,'Time_(_m_s_)')
    ax.FontSize = 30;
    configureFigureToPDF (fig);
    figName = ['Unit ', num2str(ID{a}), ' Waveform'];
    savefig(fig,fullfile(WaveformDir, [figName, '.fig']), 'compact');
end
end