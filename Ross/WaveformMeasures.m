function WvMeasures = WaveformMeasures(ID, clWaveforms, fs)
WvMeasures = struct([]);
inds = find(ismember(clWaveforms(:,1), ID));
for a = 1:length(ID)
    fig = figure('color', 'white');
    ind = inds(a);
    bartho = getBarthoMeasure((mean(clWaveforms{ind,2},2)), fs);
    plot(mean(clWaveforms{ind,2},2), 'Color', 'b');
    fig = gcf;
    ax = fig(1).Children;
    sz = size(clWaveforms{ind,2}(:,1));
    sz = sz(1);
    
    trough = min(mean(clWaveforms{ind,2},2));
    troughInd = find(ismember(mean(clWaveforms{ind,2},2), trough));
    peak2 = max(mean(clWaveforms{ind,2}(troughInd:end,:),2));
    peak2Ind = find(ismember(mean(clWaveforms{ind,2},2), peak2));
    hold on
    plot(peak2Ind*ones(100), linspace(trough, peak2, 100),...
        'LineStyle', '--', 'Color', 'k');
    plot(linspace(troughInd, peak2Ind,peak2Ind-troughInd), ...
        trough*ones(peak2Ind-troughInd), 'Color', 'k');
    text(peak2Ind + 1, trough, [num2str(bartho(1)*1e3), '_{ms}'])
    hold off
    
    
    
%     ax.XLim = [1, sz];
%     ax.YLim = [1.1*trough, -1.1*trough]; % round here
%     ax.YTick = round(linspace(ax.YLim(1), ax.YLim(2), 10));
%     ax.YTickLabel = round(ax.YTick./1e3,2);
%     ax.XTick = linspace(1,sz,10);
%     ax.XTickLabel = round(ax.XTick./fs*1e3,1);
%     ax.XTickLabelRotation = -45;
%     title(['Unit ',clWaveforms{ind,1},' Waveform']);
%     xlabel(ax,'Time_{(ms)}');
%     ylabel(ax,'Amplitude_{(mV)}');
%     ax.FontSize = 20;
    if class(ID) == 'cell'
        figureName = ID{a};
    else
        figureName = ID;
    end
    WvMeasures(a).ID = figureName;
    WvMeasures(a).PeakTrough = bartho(1);
    WvMeasures(a).HalfWidth = bartho(2);
    WvMeasures(a).FirstPeak = bartho(3);
    WvMeasures(a).SecondPeak = bartho(4);
    if class(ID) == 'char'
        return
        
    end
    
end
end