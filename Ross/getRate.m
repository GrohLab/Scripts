function rates = getRate(ID, binSize, sortedData, numSamples, fs)
Nt = numSamples/fs;
% Useless clusters (labeled as noise or they have very low firing rate)
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
bads = find(badsIdx);
totSpkCount = cellfun(@numel,sortedData(:,2));
clusterSpikeRate = totSpkCount/Nt;
silentUnits = clusterSpikeRate < 0.1;
bads = union(bads,find(silentUnits));
goods = setdiff(1:size(sortedData,1),bads);
badsIdx = badsIdx | silentUnits;
spkSubs = cellfun(@(x) round(x.*fs), sortedData(goods,2),...
    'UniformOutput', false);
gclID = sortedData(goods,1);


if binSize == round(binSize) || binSize >= 1
    binSize = binSize*(1e-3);
end

binSamples = fs*binSize;
nBins = round(numSamples/binSamples);
nUnits = length(spkSubs);
counts = zeros(nBins, nUnits);

for unitNo = 1:nUnits   % unit by unit
    unitInd = ismember(gclID,ID{unitNo});
    train = spkSubs{unitInd}'; % need to index into gclID not spkSubs
    for bin = 1:nBins   % bin by bin 
        logicalTrain = (train >(bin-1)*binSamples) & (train <= bin*binSamples);
        counts(bin,unitNo) = sum(logicalTrain);
    end
end

rates = counts/binSize;