function [burstSpkFreq, burstCF, nBursts, nSpikes, eventRatio] = getBurstingMeasures(UnitID, ISIcutoff, Traincutoff, minSpksperBurst, sortedData, samplingFrequency)

% Returns the frequency of bursts with element(n) number of spikes in (up
% to 25 spikes per burst), as well as the bursting coefficient (no of burst
% spikes relative to total spikes) and the number of bursts.
% spb = minimum spikes per burst
ind = ismember(sortedData(:,1), UnitID);
SpikeTrain = sortedData{ind,2};


fs = samplingFrequency;
rP = 10^-3;
if Traincutoff == round(Traincutoff)
    Traincutoff = Traincutoff*(10^-3);
end

if ISIcutoff == round(ISIcutoff)
    ISIcutoff = ISIcutoff*(10^-3);
end

dim=size(SpikeTrain); if dim(2)>dim(1),SpikeTrain=SpikeTrain';end  %consistent column of spike times.
if  SpikeTrain(1) == round(SpikeTrain(1))
    SpikeTrain = SpikeTrain/fs;
end
SpikeTrain(diff(SpikeTrain) < rP) = [];
ISIs = diff(SpikeTrain);
Events = find(ISIs > ISIcutoff);
Events = [1;Events];
Events=unique(Events);



c = 1;
Bursts = [];
for a = 1:length(Events) - 1
    current = Events(a); next = Events(a+1);
    if  sum(ISIs(current:next-1)) <= Traincutoff && next-current-minSpksperBurst >= 0 % min number of spikes per burst has to be greater than 1
        Bursts{c} = (SpikeTrain(current+1:next));
        c = c + 1;
    end
end
current = Events(end);

if  sum(ISIs(current:end)) <= Traincutoff && next-current-minSpksperBurst >= 0 % min number of spikes per burst has to be greater than 1
        Bursts{c} = (ISIs(current:end));
end
if length(Bursts) == 0
    nBursts = 0;
    burstSpkFreq = zeros(25,1);
    burstCF = 0;
else
    nBursts = length(Bursts);
    
    eventRatio = nBursts/(length(events)-nBursts);
    
    burstSpikes = cat(1, Bursts{:});
    burstCF = length(burstSpikes)/length(SpikeTrain);
    
    burstSz = zeros(size(Bursts))';
    for a = 1:length(Bursts)
        burstSz(a,1) = numel(Bursts{a});
    end
    burstSpkFreq = zeros(25,1);
    for a = 1: 25
        burstSpkFreq(a,1) = sum(burstSz == a);
    end
end
nSpikes = length(SpikeTrain);
end
