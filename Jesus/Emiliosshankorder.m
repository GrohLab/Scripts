dataDir
clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
clInfo
clInfo = getClusterInfo('D:\191211ProbeE1_Jesus_Emilio_Jittering_3470_1500_1500_3mW\cluster_info.tsv');
clInfo
dataDir = 'D:\191211ProbeE1_Jesus_Emilio_Jittering_3470_1500_1500_3mW';
fullfile('path','to','your','file')
chanMap = readNPY(fullfile(dataDir,'channel_map.npy'));
chanPos = readNPY(fullfile(dataDir,'channel_positions.npy'));
figure; scatter(chanPos(:,1), chanPos(:,2)); text(chanPos(:,1), chanPos(:,2), arrayfun(@num2str, chanMap,'UniformOutput',0))
[m,b] = lineariz(chanPos(:,1),4,1);
shankInd = round(m*chanPos(:,1) + b)
clInfo
chanMap
chanMap(shankInd == 1)
chShank1 = chanMap(shankInd == 1);
chShank2 = chanMap(shankInd == 2);
chShank3 = chanMap(shankInd == 3);
chShank4 = chanMap(shankInd == 4);
clInfo.shank(ismember(clInfo.channel,chShank1)) = 1
clInfo.shank(ismember(clInfo.channel,chShank2)) = 2
clInfo.shank(ismember(clInfo.channel,chShank3)) = 3;
clInfo.shank(ismember(clInfo.channel,chShank4)) = 4;
baseName = fullfile(dataDir,'Jittering_3mW_2');
save([baseName, '.mat'],'clInfo');
save([baseName, 'cluster_info.mat'],'clInfo');
save([baseName, '_clusterInfo.mat'],'clInfo');
clInfo
ismember(clInfo.channel,'7')
ismember(clInfo.channel,7)
clInfo(1,1)
clInfo(4:7,1:2)
clInfo('81','Amplitude')
clInfo.shank == 1