function addWaveformMeasures(clInfo, ID, clWaveforms, dataDir, fs)

for a = 1:length(ID)
    ind = find(ismember(clWaveforms(:,1), ID{a}));
    bartho(a,:) = getBarthoMeasure((mean(clWaveforms{ind,2},2)), fs);
  
end
clInfo = addvars(clInfo,zeros(height(clInfo),1),'Before','firing_rate',...
'NewVariableNames','Trough_Peak');
clInfo = addvars(clInfo,zeros(height(clInfo),1),'Before','firing_rate',...
'NewVariableNames','Half_Peak_Width');

inds = find(ismember(clInfo.id, ID));
clInfo.Trough_Peak(inds) = bartho(:,1);
clInfo.Half_Peak_Width(inds) = bartho(:,2);
writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'));
end