function fig = PeakTrough(clWaveforms, fs)
PT = zeros(length(clWaveforms),1);
for a = 1: length(clWaveforms)
    avWv{a} = mean(clWaveforms{a,2}')';
    
    minVal = min(avWv{a});
    minInd = find((avWv{a} == minVal));
    maxVal = max(avWv{a});
    maxInd = find((avWv{a} == maxVal));
    SampleDif = maxInd - minInd;
    msecDif = SampleDif/fs*1000;
    PT(a,1) = msecDif;
end
xLow = min(PT);
xHigh = max(PT);
scaled = round((xHigh-xLow)*20);
fig = figure('color', 'white'); histogram(PT,scaled);
title('Peak-Trough delta-t Distribution');
xlabel('Time (msec)');
end