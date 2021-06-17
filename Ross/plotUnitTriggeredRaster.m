function plotUnitTriggeredRaster(ID, clInfo, sortedData, CondTriggers, samplingFreq, Triggers, clr, rasclr)

% CondTriggers from Conditions struct



% % Redefining the stimulus signals from the low amplitude to logical values
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
% continuousNameSub(continuousNameSub == 0) = [];
% trigNames = trigNames(continuousNameSub);


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


[~, Nt, NTa] = size(discStack);
% Computing the time axis for the stack
tx = (0:Nt - 1)/fs + timeLapse(1);

% delayFlags
delayFlags = true(NTa,1);
Na = sum(delayFlags,1);
trigSubset = sort(randsample(Na,Na));
tLoc = find(delayFlags);
tSubs = tLoc(trigSubset);
% Trigger subset for stimulation shading
trigAlSubs = CondTriggers(trigSubset,:);
timeDur = round(diff(trigAlSubs, 1, 2)/fs, 3);
trigChange = find(diff(timeDur) ~= 0);



binSz = 0.01;
if class(ID) == 'cell'
    figureName = ID{1};
else
    figureName = ID;
end
% figure('Name', ['Unit ', figureName], 'Color', 'White');


ax = gca;
plotTaggedRasterFromStack(discStack([1,2],:,tSubs),...
    timeLapse, fs,figureName, ax, rasclr);
ax.YTick = [ax.YLim(1), ax.YLim(2)];
% PSTH

[~, trig, ~] = getPSTH(discStack([1,2],:,:),timeLapse,...
    ~delayFlags,binSz,fs);
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
%% Adding TTL Traces


[stims, r, csNames] = getTriggeredTTL(cst, delayFlags, trigNames, clr, Nt);


hold on
yyaxis right
for cs = r
    if exist('IDs','var')
        plot(trigTX,stims(:,cs),'LineStyle',':','LineWidth',1.5,...
            'DisplayName', IDs{cs}, 'Color', clr)
    else
        plot(trigTX,stims(:,cs),'LineStyle',':','LineWidth',1.5, 'Color', clr)
    end
    
end
% fInd = [1:length(csNames)] - 1;
f = get(gca,'Children');
lgd = legend(f((end-0)), csNames);
ax = gca;
ax.Box = 'off';

ax.YAxis(2).Visible = 'off';




ax.FontSize = 8;
ax.FontName = 'Arial';
end

