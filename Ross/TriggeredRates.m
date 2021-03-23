function Rates = TriggeredRates(Array, TriggerTimes, samplingFrequency, timeBin)

% TriggeredRates
% Feed in Unit spike times (1 cell per unit e.g. sortedData(goods,2) and trigger times (e.g. Conditions(2).Triggers)...
%...get Latencies in ms

Rates = cell(length(Array),1);
triggers = TriggerTimes/samplingFrequency;

for ind = 1:length(Array)
    Train = Array{ind};
    rate = zeros(1,length(triggers));
    
    for TrigNo = 1:length(triggers)
        
        spkInds = find(Train > triggers(TrigNo,1) & Train < triggers(TrigNo,1) + timeBin);
        
        if isempty(spkInds)
            rate(TrigNo) = 0;
        else
            cts = length(spkInds);
            rate(TrigNo) = cts/timeBin;
        end
    end
    Rates{ind} = rate;
    
end
end


