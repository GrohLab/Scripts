function plotUnitTriggeredPSTH(ID, clInfo, sortedData, CondTriggers, samplingFreq, Triggers, clr)

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
% Trigger subset for stimulation shading
trigAlSubs = CondTriggers(trigSubset,:);
timeDur = round(diff(trigAlSubs, 1, 2)/fs, 3);
trigChange = find(diff(timeDur) ~= 0);



binSz = 0.05;
if class(ID) == 'cell'
    figureName = ID{1};
else
    figureName = ID;
end
% figure('Name', ['Unit ', figureName], 'Color', 'White');

% subplot(2,1,1)
% ax = gca;
% plotRasterFromStack(discStack([1,2],:,tSubs),...
%     timeLapse, fs,figureName, ax);
% hold on
% initSub = 0;
% optsRect = {'EdgeColor',clr,'FaceColor','none'};
% for ctr = 1:numel(trigChange)
%     rectangle('Position',[0, initSub,...
%         timeDur(trigChange(ctr)), trigChange(ctr)],optsRect{:})
%     initSub = trigChange(ctr);
% end
% rectangle('Position', [0, initSub, timeDur(Na),...
%     Na - initSub],optsRect{:})


% PSTH




[PSTH, trig, sweeps] = getPSTH(discStack([1,2],:,:),timeLapse,...
    ~delayFlags,binSz,fs);
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));





PSTHax = linspace(timeLapse(1),timeLapse(2),size(PSTH,2));
area(PSTHax, PSTH, 'FaceColor', [1,0,0], 'LineStyle', 'none', 'FaceAlpha', 1);
ax = gca;

sc = 1/(Na*binSz);

cnv = ax.YLim(2)*sc;
maxYLabel = round(cnv + 5, -1);
mxYTick = ax.YLim(2)*maxYLabel/cnv;
ax.YLim(2) = mxYTick;
yscl = [ax.YLim(1):10/sc: ax.YLim(2)];
ax.YAxis.TickValues = yscl;
ax.YAxis.TickLabels = sc.*ax.YAxis.TickValues;


ax.XLim = timeLapse;
ax.YAxis.Label.String = ['Firing Rate _{(Hz)}'];
ax.FontSize = 8;
ax.FontName = 'Arial';
hold on

[stims, r, csNames] = getTriggeredTTL(cst, delayFlags, trigNames, clr, Nt);
yyaxis right
for cs = r
    if exist('IDs','var')
        plot(gca, trigTX,stims(:,cs),'LineStyle',':','LineWidth',1.5,...
            'DisplayName', IDs{cs}, 'Color', clr)
    else
        plot(gca, trigTX,stims(:,cs),'LineStyle',':','LineWidth',1.5, 'Color', clr)
    end
    
end
% fInd = [1:length(csNames)] - 1;
f = get(gca,'Children');
lgd = legend(f((end-0)), csNames);
lgd.Box = 'off';
ax = gca;
ax.Box = 'off';

ax.YAxis(2).Visible = 'off';
end




