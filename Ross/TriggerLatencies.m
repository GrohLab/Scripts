function [Latencies, Fidelities] = TriggerLatencies(Array, TriggerTimes, samplingFrequency, timeBin)

% TriggerLatencies
% Feed in Unit spike times (1 cell per unit e.g. sortedData(goods,2) and trigger times (e.g. Conditions(2).Triggers)...
%...get Latencies in ms

Latencies = cell(length(Array),1);
Fidelities = cell(length(Array),1);
triggers = TriggerTimes/samplingFrequency;

for ind = 1:length(Array)
    Train = Array{ind};
    firstSpks = zeros(1,length(triggers));
    responses = 0;
    for TrigNo = 1:length(triggers)
        
        spks = find(Train > triggers(TrigNo,1) & Train < triggers(TrigNo,1) + timeBin);
        
        if isempty(spks)
            firstSpks(TrigNo) = 0;
        else
            responses = responses + 1;
            frstVal = spks(1);
            latency = Train(frstVal) - triggers(TrigNo, 1);
            firstSpks(TrigNo) = latency;
        end
    end
    firstSpks(firstSpks == false) = [];
    Latencies{ind} = firstSpks;
    Fidelities{ind} = (responses/length(triggers))*100;
end
end


