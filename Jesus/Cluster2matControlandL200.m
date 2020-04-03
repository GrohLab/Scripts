ClustersControlOnset = nan(500,size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(2).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(2).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersControlOnset (1:N,i) = Spikes;
end

save('ClustersControlOnset.mat','ClustersControlOnset')

ClustersL200Onset = nan(500,size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1));
for i = 1:size(relativeSpkTmsStruct(1).SpikeTimes(:,1),1)
    Spikes=cell2mat(relativeSpkTmsStruct(1).SpikeTimes(i,1:end));
    N= size(Spikes,2);
    ClustersL200Onset (1:N,i) = Spikes;
end
save('ClustersL200Onset.mat','ClustersL200Onset')

nbin=350

figure
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL200Onset,'Orientation','horizontal')
title('L200')

figure
subplot(2,1,1)
hist(ClustersControlOnset,nbin)
xlabel('time(s)'),ylabel('counts')
title('Control')
subplot(2,1,2)
hist(ClustersL100Onset,nbin)
xlabel('time(s)'),ylabel('counts')
title('L100')

allControl=ClustersControlOnset(:)
allL200=ClustersL200Onset(:)

subplot(2,1,1)
histfit(allControl,nbin,'kernel')
xlabel('time(s)'),ylabel('counts')
title('Control')
subplot(2,1,2)
histfit(allL200,nbin,'kernel')
xlabel('time(s)'),ylabel('counts')
title('L200')