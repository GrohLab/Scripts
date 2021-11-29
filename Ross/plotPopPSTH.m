filterIdx = true(1,3);
filterIdx = [filterIdx; idx & wruIdx];
dimID = size(filterIdx);
nFilters = dimID(2);

csNames = fieldnames(Triggers);
Nbn = diff(timeLapse)/binSz;
popClrs = linspace(0.5, 1, nFilters);

for ccond = [3, 6, 9, 12, 15] %1:Nccond
    PSTHarray = cell(nFilters, 1);
    
    for f = 1:nFilters
        PSTHarray{f} = zeros(nnz(filterIdx(:,f)) - 1, Nbn, Nccond);
        [PSTHarray{f}(:,:,ccond), trig, sweeps] = getPSTH(discStack(filterIdx(:,f),:,:),timeLapse,...
            ~delayFlags(:,ccond),binSz,fs);
    end
    fig = figure('Name', [expName, '_', consCondNames{ccond}, '_Population activity'], 'Color', 'White');
    
    
    ax1 = subplot(4,1,1:3);
    hold(ax1, 'on')
    title([expName, '_', consCondNames{ccond}, '_Population activity_', binSz, 'binSz'])
    for f = 1:nFilters
        PSTH = PSTHarray{f}(:,:,ccond);
        [Ncl, Npt] = size(PSTH);
        psthTX = linspace(timeLapse(1),timeLapse(2),Npt);
        popPSTH = sum(PSTH,1,'omitnan')/(Ncl * sweeps);
        colr = zeros(1,3);
        colr(f) = 1;
        plot(ax1, psthTX,popPSTH,'Color',colr);  %'Color',popClrs(f)*ones(1,3)
       
    end
    legend([{'Early Responders'}, {'VPL'}, {'Late Responders'}])
    
    
    ax2 = subplot(4,1,4);
    hold on
    trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));
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
    
    
    [r,c] = size(stims);
    if r < c
        stims = stims';
    end
    
    
    plot(ax2,trigTX,stims(:,cs),'LineStyle','-','LineWidth',0.5)
    legend(consCondNames{ccond})
    linkaxes([ax1,ax2],'x')
    
    
    
    
    
    hold off
    savefig(fig, fullfile('Z:\Ross\Experiments\smrxFiles\16.12.20\12.6.21\Figures\', [fig.Name, '_', num2str(binSz), '_binSz.fig']));
end