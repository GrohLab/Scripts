%% SingleTrialPSTH

% To be run after other PSTHs


for trial = 1:length(Conditions(32).Triggers)
    Conditions(38).Triggers = Conditions(32).Triggers(trial, :);
    consideredConditions = 38;
    Nccond = length(consideredConditions);
    
    
    delayFlags = false(NTa,Nccond);
    counter2 = 1;
    for ccond = consideredConditions
        delayFlags(:,counter2) = ismember(Conditions(chCond).Triggers(:,1),...
            Conditions(ccond).Triggers(:,1));
        counter2 = counter2 + 1;
    end
    Na = sum(delayFlags,1);
    
    
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
        [PSTH(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx,:,:),timeLapse,...
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
        psthFigs(ccond) = plotClusterReactivity(PSTH(ordSubs,:,ccond),trig,sweeps,timeLapse,binSz,...
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
    
end
