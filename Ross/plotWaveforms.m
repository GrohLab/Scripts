function fig = plotWaveforms(clWaveforms)
Wv = clWaveforms(:,2);
Wv = cellfun(@transp,Wv , 'UniformOutput', false);
mnWv = cellfun(@mean,Wv, 'UniformOutput', false);
mnWv = cell2mat(mnWv);
dimWv = size(mnWv);
nCl = dimWv(1);
fig = figure;
hold on
for unit = 1:nCl
    plot(mnWv(unit,:));
end
end
