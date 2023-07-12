% function [burstSpkFreq, burstCF, nBursts, nSpikes, eventRatio] = getBurstingMeasuresMUA(UnitID, SpikeTrains, TrialStarts, TimeBefore, TimeAfter, Condition)

%%

%% Picking out the desired units
if ~exist("spkSubs")
    ind = ismember(sortedData(:,1), gclID);
    Spikes = cellfun(@(x) round(x.*fs),sortedData(ind,2),...
    'UniformOutput',false);
else
    Spikes = spkSubs;
end


%% Parameters

minSpksperBurst = 2;
Traincutoff = 40 * 10^-3;
ISIcutoff = 6 * 10^-3;
nSpkscutoff = 5;

timeBeforesecs = 2.25;
timeAftersecs = 2.25;

consconds = [20, 22, 25];


%% Triggered Spike Times


conscondnames = {Conditions(consconds).name};

Triggers = {Conditions(consconds).Triggers};
timeBefore = timeBeforesecs*fs;
timeAfter = timeAftersecs*fs;

for t = 1:length(Triggers)
    triggers = Triggers{t}(:,1);

spont_window = [triggers-timeBefore, triggers-0.01];
ev_window = [triggers+0.01, triggers+timeAfter];

consSpks = cell(length(Spikes), 2);

for cu = 1:length(Spikes)
    spks = Spikes{cu};
    sp_spks = [];
    ev_spks = [];
    for ct = 1:length(triggers)
        sp_spks = [sp_spks; spks(spks > spont_window(ct,1) & spks < spont_window(ct,2))];
        ev_spks = [ev_spks; spks(spks > ev_window(ct,1) & spks < ev_window(ct,2))];
    end
    consSpks{cu,1} = sp_spks;
    consSpks{cu,2} = ev_spks;
end


%% Spont and Evoked Bursting Measures
% burstSpkFreq = cell(size(consSpks));
burstCF = cell(size(consSpks)); % no. of spikes belonging to bursts / total spikes
nBursts = cell(size(consSpks)); % no. of bursts
nSpikes = cell(size(consSpks)); % no. of total spikes
bursteventRatio = cell(size(consSpks)); % no. of burst events




refperiod = 0.001;

for cu = 1:length(Spikes)
    for ct = 1:2
        spks = consSpks{cu,ct};
        if isempty(spks)
%             burstSpkFreq{cu,ct} = 0;
            burstCF{cu,ct} = 0;
            nBursts{cu,ct} = 0;
            nSpikes{cu,ct} = 0;
            bursteventRatio{cu,ct} = 0;
        else

        dim=size(spks); if dim(2)>dim(1),spks=spks';end  %consistent column of spike times.
        if  spks(1) == round(spks(1))
            spks = spks/fs;
        end
        spks(diff(spks) < refperiod) = [];
        ISIs = diff(spks);
        Events = [1; find(ISIs > ISIcutoff)];
        Events=unique(Events);



        c = 1;
        Bursts = [];
        for a = 1:length(Events) - 1
            current = Events(a); next = Events(a+1);
            if  sum(ISIs(current:next-1)) <= Traincutoff && next-current-minSpksperBurst >= 0 && next-current <= nSpkscutoff-1 % min number of spikes per burst has to be greater than 1
                Bursts{c} = (spks(current+1:next));
                c = c + 1;
            end
        end
        current = Events(end);

        if  sum(ISIs(current:end)) <= Traincutoff && next-current-minSpksperBurst >= 0 && next-current <= nSpkscutoff-1 % min number of spikes per burst has to be greater than 1
            Bursts{c} = (ISIs(current:end));
        end
        if isempty(Bursts)
            nBursts{cu,ct} = 0;
%             burstSpkFreq{cu,ct} = zeros(10,1);
            burstCF{cu,ct} = 0;
            bursteventRatio{cu,ct} = 0;
        else
            nBursts{cu,ct} = length(Bursts);

            bursteventRatio{cu,ct} = length(Bursts)/(length(Events));

            burstSpikes = cat(1, Bursts{:});
            burstCF{cu,ct} = length(burstSpikes)/length(spks);

            burstSz = zeros(size(Bursts))';
            for a = 1:length(Bursts)
                burstSz(a,1) = numel(Bursts{a});
            end
%             burstSpkFreq = zeros(10,1);
            for a = 1: 10
%                 burstSpkFreq(a,1) = sum(burstSz == a);
            end
        end
        nSpikes{cu,ct} = length(spks);

        end



    end
end
TriggeredBursts(t).ConditionName = conscondnames{t};
TriggeredBursts(t).burstCF = burstCF;
TriggeredBursts(t).bursteventRatio = bursteventRatio;
TriggeredBursts(t).nBursts = nBursts;
TriggeredBursts(t).nSpikes = nSpikes;

end
