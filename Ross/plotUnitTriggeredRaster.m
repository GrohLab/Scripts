function plotUnitTriggeredRaster(ID, clInfo, sortedData, CondTriggers, samplingFreq, Triggers, clr)

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



binSz = 0.01;
if class(ID) == 'cell'
    figureName = ID{1};
else
    figureName = ID;
end
% figure('Name', ['Unit ', figureName], 'Color', 'White');


ax = gca;
plotRasterFromStack(discStack([1,2],:,tSubs),...
    timeLapse, fs,figureName, ax);
hold on
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
% Plotting lines for depicting the on- and offset of the trigger
tObj = StepWaveform(trig,1/mean(diff(trigTX)));
tSubs = tObj.subTriggers;
if ~isempty(tSubs)
    for ct = 1:size(tSubs,1)
        for nf = 1:size(tSubs,2)
            line(ax,'XData',repmat(trigTX(tSubs(ct,nf)),1,2),'YData',[0.5,Na+0.5],...
                'Color',clr,'LineWidth',2)
        end
    end
% else
    trig = logical(trig/max(trig(:)));
    line(ax,'XData',trigTX,'YData',(Na+1.5)*trig - 0.5,...
        'Color',clr,'LineWidth',2)
end
ax.XLim = timeLapse;
ax.YAxis.Label.String = ['# Trials'];
ax.FontSize = 20;
ax.FontName = 'Arial';

% subplot(2,1,2)
% 
% 
% 
% PSTHax = linspace(timeLapse(1),timeLapse(2),size(PSTH,2));
% area(PSTHax, PSTH, 'FaceColor', [1,0,0], 'LineStyle', 'none', 'FaceAlpha', 0.1);
% ax = gca;
% 
% sc = 1/(Na*binSz);
% 
% cnv = ax.YLim(2)*sc;
% maxYLabel = round(cnv + 5, -1);
% mxYTick = ax.YLim(2)*maxYLabel/cnv;
% ax.YLim(2) = mxYTick;
% yscl = [ax.YLim(1):10/sc: ax.YLim(2)];
% ax.YAxis.TickValues = yscl;
% ax.YAxis.TickLabels = sc.*ax.YAxis.TickValues;
% 
% line('XData',trigTX,'YData',(0.25*mxYTick)*trig - 0.5,...
%     'Color',clr,'LineWidth',2);
% 
% ax.YAxis.Label.String = ['Firing Rate _{(Hz)}'];
% ax.XLabel.String = 'Time _{(s)}';
% ax.FontSize = 20;
% ax.FontName = 'Arial';






% stims = mean(cst(:,:,delayFlags,3));
%     stims = stims - median(stims,2);
%     for cs = 1:size(stims,1)
%         if abs(log10(var(stims(cs,:),[],2))) < 13
%             [m,b] = lineariz(stims(cs,:),1,0);
%             stims(cs,:) = m*stims(cs,:) + b;
%         else
%             stims(cs,:) = zeros(1,Nt);
%         end
%     end
% ax3 = subplot(totlX,1,5,'Parent',fig);
% [r,c] = size(stims);
% if r < c
%     stims = stims';
% end
% 
% for cs = 1:min(r,c)
%     if exist('IDs','var')
%         plot(ax3,trigTX,stims(:,cs),'LineStyle','-.','LineWidth',0.5,...
%             'DisplayName', IDs{cs})
%     else
%         plot(ax3,trigTX,stims(:,cs),'LineStyle','-.','LineWidth',0.5)
%     end
%     
%     ax3.Children(1).Color = defineColorForStimuli(IDs(cs));
%     
%     if cs == 1
%         ax3.NextPlot = 'add';
%     end
% end
% legend(ax3,'show','Location','best')
% 
% ax3.Box = 'off';
% ax3.XLim = [timeLapse(1), timeLapse(2)];
% ax3.XAxis.Visible = 'off';
% ax3.YAxis.Visible = 'off';
% linkaxes([ax1,ax2,ax3],'x')
% end


