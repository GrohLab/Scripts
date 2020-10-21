% PrepTables
clInfo = getClusterInfo(fullfile(dataDir,'cluster_info.tsv'));
sZ = height(clInfo);

for a  = 1:sZ
clInfo.Properties.RowNames{a,1} = [expName, '_', clInfo.id{a}];
clInfo.id{a,1} = clInfo.Properties.RowNames{a};

end
Region = cell(sZ,1);
Model = cell(sZ,1);
for a = 1:sZ
    Region{a} = 'VPL';
    Model{a} = 'CFA';
end
clInfo = addvars(clInfo,Model,'NewVariableNames','Model','After','id');
clInfo = addvars(clInfo,Region,'NewVariableNames','Region','After','id');
% writeClusterInfo(clInfo, fullfile(dataDir,[expName, '_Table'], '.tsv'));
