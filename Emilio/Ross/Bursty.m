function [BurstStruct, popFreq, totalBursts, totalSpikes] = Bursty(SpikeTrainArray, fs)
for a = 1:length(SpikeTrainArray)
    [burstSpkFreq, burstCF, nBursts, nSpikes] = getBurstingMeasures(SpikeTrainArray{a, 1}, fs, 4, 100);
    BurstStruct(a).Freq = burstSpkFreq;
    BurstStruct(a).Coeff = burstCF;
    BurstStruct(a).nBursts = nBursts;
    BurstStruct(a).nSpikes = nSpikes;
end

freq = cat(2, BurstStruct(:).Freq);
popFreq = sum(freq')';
bursts = cat(2, BurstStruct(:).nBursts);
totalBursts = sum(bursts')';
relFreq = popFreq/totalBursts;
spikes = cat(2, BurstStruct(:).nSpikes);
totalSpikes = sum(spikes')';






end