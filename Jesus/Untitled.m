ClustersControlOnsetMOCK = nan(2000,size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(1).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnsetMOCK (1:N,i) = Spikes;
end


ClustersL200OnsetMOCK = nan(2000,size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(2).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL200OnsetMOCK (1:N,i) = Spikes;
end