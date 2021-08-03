function plotUnitTriggeredResponse(ID, clInfo, sortedData, CondTriggers, samplingFreq, Triggers)

% CondTriggers from Conditions struct



% Redefining the stimulus signals from the low amplitude to logical values
whStim = {'piezo','whisker','mech','audio'};
cxStim = {'laser','light'};
lfpRec = {'lfp','s1','cortex','s1lfp'};
trigNames = fieldnames(Triggers);
numTrigNames = numel(trigNames);
ctn = 1;
continuousSignals = cell(numTrigNames,1);
continuousNameSub = zeros(size(trigNames));
while ctn <= numTrigNames
    if contains(trigNames{ctn},whStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},cxStim,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    if contains(trigNames{ctn},lfpRec,'IgnoreCase',true)
        continuousSignals{ctn} = Triggers.(trigNames{ctn});
        continuousNameSub(ctn) = ctn;
    end
    ctn = ctn + 1;
end
continuousSignals(continuousNameSub == 0) = [];
continuousNameSub(continuousNameSub == 0) = [];
trigNames = trigNames(continuousNameSub);


% Preparing stack inputs
fs = samplingFreq;
Ns = min(structfun(@numel,Triggers));
timeLapse = [-2, 7];
clInd = ismember(clInfo.id, ID);
spkLog = StepWaveform.subs2idx(round(sortedData{clInd,2}*fs),Ns);
spks = sortedData(clInd,2);

spkSubs = cellfun(@(x) round(x.*fs), spks,...
    'UniformOutput', false);

% Making discStack
[discStack, cst] = getStacks(spkLog,CondTriggers,'on',...
    timeLapse,fs,fs,spkSubs,continuousSignals);


[Ne, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);

% delayFlags
delayFlags = true(NTa,1);
Na = sum(delayFlags,1);
trigSubset = sort(randsample(Na,Na));
tLoc = find(delayFlags);
tSubs = tLoc(trigSubset);
binSz = 0.1;


plotRasterFromStack(discStack([1,2],:,tSubs),...
    timeLapse, fs,'');
ax1 = gca;
hold on
clr = [0,1,1];

[PSTH, trig, sweeps] = getPSTH(discStack([1,2],:,:),timeLapse,...
    ~delayFlags,binSz,fs);
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));

% Plotting lines for depicting the on- and offset of the trigger
tObj = StepWaveform(trig,1/mean(diff(trigTX)));
tSubs = tObj.subTriggers;
if ~isempty(tSubs)
    for ct = 1:size(tSubs,1)
        for nf = 1:size(tSubs,2)
            line(ax1,'XData',repmat(trigTX(tSubs(ct,nf)),1,2),'YData',[0.5,Na+0.5],...
                'Color',clr,'LineWidth',2)
        end
    end
else
    trig = logical(trig/max(trig(:)));
    line(ax1,'XData',trigTX,'YData',(Na+1.5)*trig - 0.5,...
       'Color',clr,'LineWidth',2)
end

yyaxis('right')
sc = 1/(Na*binSz);


PSTHax = linspace(timeLapse(1),timeLapse(2),size(PSTH,2));
area(PSTHax, PSTH, 'FaceColor', [0,0,0], 'LineStyle', 'none', 'FaceAlpha', 0.1);
ax = gca;
adjustFactor = round(sc.*ax.YAxis(2).TickValues)/sc;
ax.YAxis(2).TickValues = adjustFactor;
ax.YAxis(2).TickLabels = sc.*ax.YAxis(2).TickValues;
end


