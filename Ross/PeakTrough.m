%function [PT, fig] = PeakTrough(clWaveforms, ID, fs)
PT = zeros(length(ID),1);
r = 1;
for a = 1:length(ID)
    ind = find(ismember(clWaveforms(:,1), ID(a)));
    avWv{r} = mean(clWaveforms{ind,2}')';
    minVal = min(avWv{r});
    minInd = find((avWv{r} == minVal));
    maxVal = max(avWv{r}(minInd:end,1));
    maxInd = find((avWv{r} == maxVal));
    SampleDif = maxInd - minInd;
    msecDif = SampleDif/fs*1000;
    PT(r,1) = msecDif;
    r = r + 1;
end
xLow = min(PT);
xHigh = max(PT);
scaled = round((xHigh-xLow)*20);
fig = figure('color', 'white'); histogram(PT,scaled);
title('Peak-Trough delta-t Distribution');
xlabel('Time (msec)');
fig = gcf;
ax = fig(1).Children;
ax.FontSize = 20;
%end