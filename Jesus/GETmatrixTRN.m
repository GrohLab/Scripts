ClustersL1OnsetTRN = nan(2000,size(PopSpikeTimes1,1));
for i = 1:size(PopSpikeTimes1(:,1),1)
    Spikes=cell2mat(PopSpikeTimes1(i,1:end));
    N= size(Spikes,2);
    ClustersL1OnsetTRN (1:N,i) = Spikes;
end

ClustersL10OnsetTRN = nan(2000,size(PopSpikeTimes10,1));
for i = 1:size(PopSpikeTimes10(:,1),1)
    Spikes=cell2mat(PopSpikeTimes10(i,1:end));
    N= size(Spikes,2);
    ClustersL10OnsetTRN (1:N,i) = Spikes;
end

ClustersL50OnsetTRN= nan(2000,size(PopSpikeTimes50,1));
for i = 1:size(PopSpikeTimes50(:,1),1)
    Spikes=cell2mat(PopSpikeTimes50(i,1:end));
    N= size(Spikes,2);
    ClustersL50OnsetTRN (1:N,i) = Spikes;
end

ClustersL100OnsetTRN = nan(2000,size(PopSpikeTimes100,1));
for i = 1:size(PopSpikeTimes100(:,1),1)
    Spikes=cell2mat(PopSpikeTimes100(i,1:end));
    N= size(Spikes,2);
    ClustersL100OnsetTRN (1:N,i) = Spikes;
end

ClustersL200OnsetTRN = nan(2000,size(PopSpikeTimes200,1));
for i = 1:size(PopSpikeTimes200(:,1),1)
    Spikes=cell2mat(PopSpikeTimes200(i,1:end));
    N= size(Spikes,2);
    ClustersL200OnsetTRN (1:N,i) = Spikes;
end

ClustersControlOnsetTRN = nan(2000,size(PopSpikeTimesControl,1));
for i = 1:size(PopSpikeTimesControl(:,1),1)
    Spikes=cell2mat(PopSpikeTimesControl(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnsetTRN (1:N,i) = Spikes;
end