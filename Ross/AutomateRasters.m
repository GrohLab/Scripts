% DE_Jittering needs to be unfiltered for significance for this to work!


rasterDir = fullfile(figureDir,'Rasters\');
if ~mkdir(rasterDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end
for pwr = [2, 5, 10, 15]
MchTblInd = ['Mech_Control_', num2str(pwr), 'mW_MR'];
LasTblInd = ['Laser_Control_', num2str(pwr), 'mW_LR'];
MchCondControl = ['Mech_Control_', num2str(pwr), 'mW'];
LasCondControl = ['Laser_Control_', num2str(pwr), 'mW'];
MchLasCond = ['Mech_Laser_', num2str(pwr), 'mW'];
EffectTblInd = ['Mech_Control_', num2str(pwr), 'mW_vs_Mech_Laser_', num2str(pwr), 'mW_Evoked_Response'];
TblInd = find(clInfo.(MchTblInd)); % ATM this only makes rasters that show sig control mech response
clIDind = clInfo.id(TblInd);
lngth = length(clIDind);
for a = 1:lngth
    cl = clIDind(a);
    clSel = find(ismember(pclID, cl));
    if chCond == 1
        rasCondSel = find(ismember(consCondNames, MchCondControl) | ismember(consCondNames, MchLasCond));
        label = 'Mech';
    else
        rasCondSel = find(ismember(consCondNames, LasCondControl) | ismember(consCondNames, MchLasCond));
        label = 'Laser';
    end
    rasCond = consideredConditions(rasCondSel);
    rasCondNames = consCondNames(rasCondSel);
    Nrcl = numel(clSel);
    % Reorganize the rasters in the required order.
    clSub = find(ismember(gclID, pclID(clSel)))+1;
    [rasIdx, rasOrd] = ismember(pclID(ordSubs), pclID(clSel));
    clSub = clSub(rasOrd(rasIdx));
    clSel = clSel(rasOrd(rasOrd ~= 0));
    Nma = min(Na(rasCondSel));
    rasFig = figure;
    Nrcond = length(rasCond);
    ax = gobjects(Nrcond*Nrcl,1);
    for cc = 1:length(rasCond)
        % Equalize trial number
        trigSubset = sort(randsample(Na(rasCondSel(cc)),Nma));
        tLoc = find(delayFlags(:,rasCondSel(cc)));
        tSubs = tLoc(trigSubset);
        % Trigger subset for stimulation shading
        trigAlSubs = Conditions(rasCond(cc)).Triggers(trigSubset,:);
        timeDur = round(diff(trigAlSubs, 1, 2)/fs, 3);
        trigChange = find(diff(timeDur) ~= 0);
        for ccl = 1:Nrcl
            lidx = ccl + (cc - 1) * Nrcl;
            ax(lidx) = subplot(Nrcond, Nrcl, lidx);
            title(ax(lidx),sprintf('%s cl:%s',rasCondNames{cc},pclID{clSel(ccl)}))
            plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                timeLapse, fs,'',ax(lidx));
            ax(lidx).YAxisLocation = 'origin';ax(lidx).YAxis.TickValues = Nma;
            ax(lidx).YAxis.Label.String = Nma;
            ax(lidx).YAxis.Label.Position =...
                [timeLapse(1)-timeLapse(1)*0.65, Nma,0];
            ax(lidx).XAxis.TickLabels =...
                cellfun(@(x) str2num(x)*1e3, ax(lidx).XAxis.TickLabels,...
                'UniformOutput', 0);
            xlabel(ax(lidx), 'Time [ms]')
            initSub = 0;
            optsRect = {'EdgeColor','none','FaceColor','none'};
            for ctr = 1:numel(trigChange)
                rectangle('Position',[0, initSub,...
                    timeDur(trigChange(ctr)), trigChange(ctr)],optsRect{:})
                initSub = trigChange(ctr);
            end
            rectangle('Position', [0, initSub, timeDur(Nma),...
                Nma - initSub],optsRect{:})
        end
    end
    linkaxes(ax,'x')
    rasFigName = ['Unit_', cell2mat(cl), '_', label, '_', num2str(pwr), 'mW_Raster']; 
    configureFigureToPDF (rasFig);
    savefig(rasFig,fullfile(rasterDir, [rasFigName, '.fig']));
end
end
close all