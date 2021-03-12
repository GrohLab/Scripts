function  iOK = KobayashiPrep(dataDir, sortedData, Conditions, Triggers, fs)
iOK = false;


% Creating the Kobayashi directory
KobaDir = fullfile(dataDir,'Kobayashi\');
if ~mkdir(KobaDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end

%% Finding goods

% Number of total samples
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;
% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);

%% Finding Spontaneous activity of good clusters
spkSubs = cellfun(@(x) round(x.*fs), sortedData(goods,2),...
'UniformOutput', false);
[FR, ~, sponSpks, ~] = getSpontFireFreq(spkSubs, Conditions(1).Triggers, [0, Nt], fs, 2);
%% Saving Spikes in Kobayashi format
sponSpks = cellfun(@(x) (x/fs)*1e3, sponSpks,...
'UniformOutput', false);
formatSpec = '%f \n';
ind = FR > 0.5;

sponSpks = sponSpks(find(ind));
for a = 1:length(sponSpks)
fileID = fopen([KobaDir, ['cell', num2str(a-1)], '.txt'], 'w');
fprintf(fileID, formatSpec, sponSpks{a});
fclose(fileID);
end

iOK = true;
end
