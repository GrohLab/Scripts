%% Rasters
% DE_Jittering needs to be unfiltered for significance for this to work!
csNames = fieldnames(Triggers);
% csNames = csNames(2:end);
IDs = csNames;
trigTX = linspace(timeLapse(1),timeLapse(2),size(trig,2));

rasterDir = fullfile(figureDir,'Rasters\');
if ~mkdir(rasterDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end


Power = NaN(length(consCondNames),1);
for cc = 1:length(consCondNames)
    
    mWfind = strfind(consCondNames{cc}, 'mW');
    
    pwr = consCondNames{cc}(mWfind-2:mWfind+1);
    if isnan(pwr)
        pwrMissing = true;
    elseif contains(pwr, '.')
        pwr = consCondNames{cc}(mWfind-3:mWfind+1);
    elseif contains(pwr, '_') || contains(pwr, ' ')
        pwr = pwr(2:end);
    end
    pwr = str2double(pwr(1:end-2));
    Power(cc) = pwr;
end
pwrs = unique(Power);

med = 0.2;
        low = 0.4;
        high = 0;
colours = ones(1,3);
colours(1:4,:) = med; colours(5:8,:) = low; colours(9:10,:) = high;

for a = 1%:length(pwrs)
    pwr = pwrs(a);
    % MchTblInd = ['Mech_Control_', num2str(pwr), 'mW_MR'];
    % LasTblInd =  [20,21,22];%['Laser_Control_', num2str(pwr), 'mW_LR'];
    % MchCondControl = ['Mech_Control_', num2str(pwr), 'mW'];
    % LasCondControl = ['Laser_Control_', num2str(pwr), 'mW'];
    % MchLasCond = ['Mech_Laser_', num2str(pwr), 'mW'];
    % EffectTblInd = ['Mech_Control_', num2str(pwr), 'mW_vs_Mech_Laser_', num2str(pwr), 'mW_Evoked_Response'];
    % TblInd = find(clInfo.ActiveUnit); % ATM this only makes rasters that show sig control mech response
    % clIDind = clInfo.id(TblInd);
    pwrInd = Power == pwr;
    clIDind =  ind;
    lngth = length(clIDind);
    for a = 1:lngth
        rng('default');
        cl = clIDind(a);
        clSel = find(ismember(pclID, cl));
        %     if chCond == 1
        %         rasCondSel = find(ismember(consCondNames, MchCondControl) | ismember(consCondNames, MchLasCond));
        %         label = 'Mech';
        %     else
        %         rasCondSel = find(ismember(consCondNames, LasCondControl) | ismember(consCondNames, MchLasCond));
        %         label = 'Laser';
        %     end
        rasCondSel = [1 2 3];
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
        columns = length(pwrs);
        Nrcond = length(rasCond);
        ax = gobjects(4*Nrcond*Nrcl,1);
        lidx = 1;
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
                
                stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
                %                 stims = stims([2,3],:);
                
                stims = stims - median(stims,2);
                
                
                
                for cs = 2 %for cs = 1:size(stims,1)
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
                    stmClr = zeros(r, 3);
                    
                end
                
                
                
                
                for cs = 2 %for cs = 1:size(stims,1)
                    if stims(1,cs) > 0.5
                        stim = ones(size(stims(:,cs)))-stims(:,cs);
                    else
                        stim = stims(:,cs);
                    end
%                     dP=[0; diff(stim)];
                    stim = smooth(stim,10);
                    ax(lidx) = subplot(4*Nrcond, Nrcl, lidx);
                    if exist('IDs','var')
                        plot(trigTX,stim, 'LineStyle','-','LineWidth', 1,...
                            'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
%                         P = stim'; %ax(lidx).Children.YData;
%                         dP = [0 diff(P)];
%                         dP = smooth(dP,100);
%                         dP = dP-mean(dP(1:100));dP=dP/max(dP);
%                         
%                         yyaxis right
%                         plot(trigTX, dP, 'LineStyle','-','LineWidth', 1,...
%                             'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
                    else
                        plot(trigTX,stim,'LineStyle','-','LineWidth',1, 'Color', stmClr(cs,:))
                        
%                         P = ax(lidx).Children.YData;
%                         dP = [0 diff(P)];
%                         dP = smooth(dP,100);
%                         dP = dP-mean(dP(1:100));dP=dP/max(dP);
%                         yyaxis right
%                         plot(trigTX, dP, 'LineStyle','-','LineWidth', 1,...
%                             'DisplayName', IDs{cs}, 'Color', stmClr(cs,:))
                    end
                    ax(lidx).Visible = 'off';
                    
                    
                    %                                        ax2.Children(1).Color = defineColorForStimuli(IDs(cs));
                    
                    if cs == 1
                        ax2.NextPlot = 'add';
                    end
                    
                end
                %                 ax = gca;
                %
                %
                % %                 ax.YAxis(2).Limits = [0.015, 1];
                % %                 ax.YAxis(2).Visible = 'off';
                %                 ax.FontName ='Arial';
                %                 ax.FontSize = 12;
                
                %                 f=get(gca,'Children');
                %                 legend(f)
                %
                lidx = lidx + 1;
                
                
                %                 lidx = ccl + (cc - 1) * Nrcl;
                ax(lidx) = subplot(4*Nrcond, Nrcl, lidx:lidx+2);       %  subplot( Nrcl, Nrcond, lidx); % to plot the other way around
                title(ax(lidx),sprintf(rasCondNames{cc}), 'Interpreter', 'none') % ,pclID{clSel(ccl)}
                plotRasterFromStack(discStack([1,clSub(ccl)],:,tSubs),...
                    timeLapse, fs,'',colours(lidx,:),ax(lidx));
                ax(lidx).YAxisLocation = 'origin';ax(lidx).YAxis.TickValues = Nma;
                ax(lidx).YAxis.Label.String = 'Trials';
                %                 ax(lidx).YAxis.Label.Position =...
                %                     [timeLapse(1)-timeLapse(1)*0.65, Nma,0];
                %             ax(lidx).XAxis.TickLabels =...
                %                 cellfun(@(x) (x)*1e3, ax(lidx).XAxis.TickValues,...
                %                 'UniformOutput', 0);
                xlabel(ax(lidx), 'Time [s]')
                initSub = 0;
                optsRect = {'EdgeColor','none','FaceColor','none'};
                for ctr = 1:numel(trigChange)
                    rectangle('Position',[0, initSub,...
                        timeDur(trigChange(ctr)), trigChange(ctr)],optsRect{:})
                    initSub = trigChange(ctr);
                end
                rectangle('Position', [0, initSub, timeDur(Nma),...
                    Nma - initSub],optsRect{:})
                
                ax(lidx).XAxis.Visible = 'off';
                ax(lidx).YAxis.Visible = 'off';
                stims = mean(cst(:,:,delayFlags(:,rasCondSel(cc))),3);
                %                 stims = stims([2,3],:);
                
                stims = stims - median(stims,2);
                
                
                
                for cs = 2 %for cs = 1:size(stims,1)
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
                    stmClr = zeros(r, 3);
                    
                end
                
                %                 stmClr =[ 1 0 0 0.25; 0 1 1 0.25];
                
                
                
                lidx = lidx + 3;
                
                
            end
        end
        
        
        
        
        
        rasConds = rasCondNames{1};
        if length(rasCondNames) > 1
            for r = 2:length(rasCondNames)
                rasConds = [rasConds, '+', rasCondNames{r}];
            end
        end
        
        
        ax(2).Title.Color = colours(2,:);
        ax(6).Title.Color = colours(6,:);
        ax(10).Title.Color = colours(10,:);
        ax(1).Children.Color = colours(1,:);
        ax(5).Children.Color = colours(5,:);
        ax(9).Children.Color = colours(9,:);
        
        
        linkaxes(ax,'x')
        rasFigName = ['Unit_', cell2mat(cl), '_', ];
        rasFig.Name = [rasFigName, '_', num2str(pwr), 'mW'];
        configureFigureToPDF (rasFig);
        set(rasFig, 'Position', get(0, 'ScreenSize')/2);
        saveas(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds,'_', num2str(timeLapse(1)), '_to_', num2str(timeLapse(2)),'.emf']));
        %savefig(rasFig,fullfile(rasterDir, [rasFigName, ' ', num2str(pwr), 'mW.fig']));
        savefig(rasFig,fullfile(rasterDir, [rasFigName,'_',rasConds, '.fig']));
    end
end

