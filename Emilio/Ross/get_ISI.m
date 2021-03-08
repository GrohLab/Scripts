%% ISI PDF & CDF
modFlags = clInfo.Mech_Control_10mW_MR; % Whatever you wanna see buddy
binFig = figure('Visible','off');
areaFig = figure('Visible','on', 'Color', [1,1,1],'Name','ISI probability');
areaAx = gobjects(2,1);
lax = log(1/fs):0.2:log(0.03);
% cmap = [0,0,102;... azul marino
%     153, 204, 255]/255; % azul cielo
cmap = lines(Nccond);
%areaOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.3, 'FaceColor'};
areaOpts = {'Color','LineStyle','--','LineWidth',1.5};
ley = ["Potentiated", "Depressed"];
isiPdf = zeros(length(lax)-1,2*Nccond);
cmod = 1; % modulated
Nccl = sum(modFlags(:,cmod));
for ccond = 1:Nccond
    cnt = [cmod-1,ccond] * [2,1]';
    ISI = cellfun(@(x) diff(x(x >= -0.035 & x <= 0)), ...
        relativeSpkTmsStruct(ccond).SpikeTimes(modFlags(:,cmod),:),...
        'UniformOutput', 0);
    ISI_merge = [ISI{:}];
    lISI = log(ISI_merge);
    hisi = histogram(lISI, lax, 'Parent', binFig, 'DisplayStyle', 'stairs');
    spkPerT_C = hisi.Values./(NaStack(ccond)*Nccl);
    isiPdf(:,cnt) = spkPerT_C/sum(spkPerT_C);
    plot(areaAx(cmod), 1e3*exp(lax(1:end-1)), isiPdf(:,cnt),...
        areaOpts{1}, cmap(ccond,:), areaOpts{4:5});
    if ccond == 1
        hold(binFig.Children, 'on'); hold(areaAx(cmod), 'on');
        areaAx(cmod).XAxis.Scale = 'log';
        areaAx(cmod).XLabel.String = "ISI [ms]";
        areaAx(cmod).YLabel.String = "ISI probability";
    end
    yyaxis(areaAx(cmod) ,'right')
    
    plot(areaAx(cmod), 1e3*exp(lax(1:end-1)), cumsum(isiPdf(:,cnt)),...
        areaOpts{1}, cmap(ccond,:), areaOpts{2:3})
    ylim(areaAx(cmod), [0,1]); areaAx(cmod).YAxis(2).Color = [0,0,0];
    ylabel(areaAx(cmod), 'Cumulative probability');
    yyaxis(areaAx(cmod) ,'left'); set(areaAx(cmod),'Box','off')
end
title(areaAx(cmod), sprintf('ISI PDF for %s clusters',ley(cmod)));
legend(areaAx(cmod),consCondNames)
hisi = hisi.Parent.Children; linkaxes(areaAx,'xy')
