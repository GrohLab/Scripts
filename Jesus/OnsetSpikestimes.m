

    ControlOnsetFlags = ClustersControlOnset >=0.0 & ClustersControlOnset <=0.05;
    PopControlOnsetTimes = ClustersControlOnset (ControlOnsetFlags);
 
OnsetSpikesTimesControl= nan(500,(max(sum(ControlOnsetFlags,1))));
for Sz=1:size(ClustersControlOnset,2)
      OnsetTemporalSpikesTimesControl=ClustersControlOnset(ControlOnsetFlags(:,1:end));
      N = size(OnsetTemporalSpikesTimesControl,1);
      OnsetSpikesTimesControl(1:N,1:Sz)= OnsetTemporalSpikesTimesControl;
      
end

ClustersControlOnset = nan(500,size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(2).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnset (1:N,i) = Spikes;
end