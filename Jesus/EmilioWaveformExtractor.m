mean_wf_cell = cellfun(@(x) mean(x,2), clWaveforms(:,2), 'UniformOutput', 0);
mean_wf = cell2mat(mean_wf_cell');size(mean_wf)
txwf = (0:size(mean_wf,1)-1)/fs;
figure; plot(txwf, mean_wf)
feats = getWaveformFeatures(mean_wf, fs);
figure; scatter(feats(:,1), feats(:,2))
text(feats(:,1), feats(:,2), clWaveforms(:,1))

wideEntropyIdx=feats(:,2)>0.5 & feats(:,2)<0.725 ;
% CUIDADO necesitas tu bandera para asignar las espigas que necesites
wideSpikeID = clWaveforms(wideEntropyIdx,1);
wruIDname = gclID(wruIdx); %conocer el nombre de las unidades responsivas, wruIdx es logico
colFiltIdx = ismember(wruIDname, wideSpikeID);
ControlWEntropyOnset=ClustersControlOnset(:,colFiltIdx)
% [wruID(colFiltIdx), wideSpikeID]

figure; scatter(feats(wideEntropyIdx,1),feats(wideEntropyIdx,2))
text(feats(wideEntropyIdx,1), feats(wideEntropyIdx,2), clWaveforms(wideEntropyIdx,1))