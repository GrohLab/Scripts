%ClInfo Depositing given that you have run DE_Jittering
Triggers = struct('MechanicalTTL',continuousSignals{1});
Hcno = Results(1).Activity(1).Pvalues < 0.05;
Hcno(:,2) = Results(1).Activity(2).Pvalues < 0.05;
clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
chanMap = readNPY(fullfile(dataDir,'channel_map.npy'));
chanPos = readNPY(fullfile(dataDir,'channel_positions.npy'));
[m,b] = lineariz(chanPos(:,1), 6, 1)
shank = m*chanPos(:,1) + b
shank = round(shank)
shMap = containers.Map(chanMap, shank);
setShank = @(x) shMap(x);
clInfo.shank = arrayfun(setShank, clInfo.channel);
groupsummary(clInfo,'group')
size(badsIdx)
clInfo = addvars(clInfo,ActiveUnit,'NewVariableNames','ActiveUnit','After','id');
clInfo.ActiveUnit(~badsIdx','NewVariableNames','ActiveUnit','After','id');
groupsummary(clInfo,'shank')
groupsummary(clInfo,'shank','nnz','ActiveUnit')
clInfo{gclID(H(:,1)),'ControlMR'} = true;
clInfo{gclID(H(:,2)),'CNO_MR'} = true;
clInfo{gclID(Hcno(:,1)),'Spont_CNO_R'} = true;
clInfo{gclID(Hcno(:,2)),'Evoked_CNO_R'} = true;
writeClusterInfo(clInfo, fullfile(dataDir,'cluster_info.tsv'))