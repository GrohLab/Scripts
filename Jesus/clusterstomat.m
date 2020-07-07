ClustersControlOnset = nan(1200,size(relativeSpkTmsStruct(6).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(6).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(6).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnset (1:N,i) = Spikes;
end

save('ClustersControlOnset.mat','ClustersControlOnset')

ClustersL200Onset = nan(1200,size(relativeSpkTmsStruct(5).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(5).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(5).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL200Onset (1:N,i) = Spikes;
end
save('ClustersL200Onset.mat','ClustersL200Onset')

ClustersL100Onset = nan(1200,size(relativeSpkTmsStruct(4).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(4).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(4).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL100Onset (1:N,i) = Spikes;
end
save('ClustersL100Onset.mat','ClustersL100Onset')

ClustersL50Onset = nan(1200,size(relativeSpkTmsStruct(3).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(3).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(3).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL50Onset (1:N,i) = Spikes;
end
save('ClustersL50Onset.mat','ClustersL50Onset')

ClustersL10Onset = nan(1200,size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(2).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL10Onset (1:N,i) = Spikes;
end
save('ClustersL10Onset.mat','ClustersL10Onset')

ClustersL1Onset = nan(1200,size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(1).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL1Onset (1:N,i) = Spikes;
end
save('ClustersL1Onset.mat','ClustersL1Onset')

