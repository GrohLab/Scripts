%% prepPopPSTH

filtStr = 'unfiltered';
%% Filtering

dimID = size(idx);
nFilters = dimID(2);
t = true(1,nFilters);
filterIdx = [t; idx];
filtStr = 'filtered';


ordSubs = 1:nnz(filterIdx(2:Ncl+1,1));
pclID = gclID(filterIdx(2:Ncl+1,1));
orderedStr = [];

%% Plot PSTH
goodsIdx = logical(clInfo.ActiveUnit);
csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
if (Nbn - round(Nbn)) ~= 0
    Nbn = ceil(Nbn);
end
PSTH = zeros(nnz(filterIdx) - 1, Nbn, Nccond);
psthFigs = gobjects(Nccond,1);
for ccond = 1:Nccond
    figFileName = sprintf('%s %s VW%.1f-%.1f ms B%.1f ms RW%.1f-%.1f ms SW%.1f-%.1f ms %sset %s (%s)',...
        expName, Conditions(consideredConditions(ccond)).name, timeLapse*1e3,...
        binSz*1e3, responseWindow*1e3, spontaneousWindow*1e3, onOffStr,...
        orderedStr, filtStr);
    [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx(:,1),:,:),timeLapse,...
        ~delayFlags(:,ccond),binSz,fs);
    stims = mean(cst(:,:,delayFlags(:,ccond)),3);
    stims = stims - median(stims,2);
    for cs = 1:size(stims,1)
        if abs(log10(var(stims(cs,:),[],2))) < 13
            [m,b] = lineariz(stims(cs,:),1,0);
            stims(cs,:) = m*stims(cs,:) + b;
        else
            stims(cs,:) = zeros(1,Nt);
        end
    end
    psthFigs(ccond) = plotPopPSTH(PSTH(ordSubs,:,ccond),trig,sweeps,timeLapse,binSz,...
        [{Conditions(consideredConditions(ccond)).name};...
        pclID(ordSubs)],...
        strrep(expName,'_','\_'),...
        stims, csNames);
    configureFigureToPDF(psthFigs(ccond));
    psthFigs(ccond).Children(end).YLabel.String =...
        [psthFigs(ccond).Children(end).YLabel.String,...
        sprintf('^{%s}',orderedStr)];
%     if ~exist([figFileName,'.pdf'], 'file')
%         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.pdf']),...
%             '-dpdf','-fillpage')
%     end
%     if ~exist([figFileName,'.emf'], 'file')
%         print(psthFigs(ccond), fullfile(figureDir,[figFileName, '.emf']),...
%             '-dmeta')
%     end
    
end
% for a = 1:length(consideredConditions)
%     savefig(figure(a), fullfile(figureDir, [consCondNames{a}, '_filtered_PSTH_0.001binSz.fig']));
% end

