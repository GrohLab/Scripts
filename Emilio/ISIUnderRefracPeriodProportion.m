fnOpts = {'UniformOutput', false};
rpTh = [0.5, 1, 2, 5, 10];

dataDir = 'Z:\Emilio\SuperiorColliculusExperiments\Anaesthetised\M24_220224_SC-E1_IO&Freq';
load(fullfile(dataDir, "SC_CVG_all_channels.mat"), "sortedData", "fs")
clInfo = readClusterInfo(fullfile(dataDir, "cluster_info.tsv"));
clIsi = cellfun(@diff, sortedData(clInfo.ActiveUnit==1,2), fnOpts{:});

isiPC = cellfun(@(x) 100*sum(x<(rpTh*1e-3))/numel(x), clIsi, fnOpts{:});
isiPC = cat(1,isiPC{:});

gclID = clInfo.Properties.RowNames(clInfo.ActiveUnit==1);
gclID = cellfun(@(x) cat(2, '5.',x), gclID, fnOpts{:});

save(fullfile(dataDir, "Refractory periods 5.mat"), "clIsi", "isiPC", ...
    "gclID", "rpTh")