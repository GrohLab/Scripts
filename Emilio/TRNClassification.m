%% TRN clustering for response
Ngcl = sum(wruIdx);
chMap = readNPY(fullfile(dataDir,'channel_map.npy'));
chPos = readNPY(fullfile(dataDir,'channel_positions.npy'));
orthMat = [0,1;-1,0];
ptsPer = 0.2;
if sum(wruIdx) < 10
    ptsPer = 0.9;
end

featMat = [-1./mdls(:,2), qDiff, -1./cmdls(wruIdx,2), cqDiff(wruIdx)].*1e3;

rmdl = boot_fit_poly(featMat(:,[1,2]), 1, ptsPer, 100);
rn = getHesseLineForm(rmdl);
rI = featMat(:,[1,2]) * orthMat * rn;

amdl = boot_fit_poly(featMat(:,[3,4]), 1, ptsPer, 100);
an = getHesseLineForm(amdl);
aI = featMat(:,[3,4]) * orthMat * an;

trnMat = [rI, aI]; [trnMat_w, whMat, trnMu] = whitenPoints(trnMat);
try
    gm = fitgmdist(trnMat_w, 2);
    trnFlag = gm.cluster(trnMat_w);
catch
    fprintf(1,'Time to classify with your eyes\n')
    figure; scatter(trnMat(:,1), trnMat(:,2), '.');
    text(trnMat(:,1), trnMat(:,2), pclID)
    return
end

trnSubs = find(ismember(gclID, pclID(trnFlag))); sqrSb = sqrt(length(trnSubs));
trnSubs_psth = find(trnFlag);
apeFig = figure('Name', 'Auto-correlation & PSTH', 'Color',[1,1,1]);
for csp = 1:size(trnSubs,1)
    subplot(ceil(sqrSb),ceil(sqrSb),csp);
    plot(corrTx, acorrs(trnSubs(csp),:))
    yyaxis right; plot(btx, PSTH(trnSubs_psth(csp),:,1))
    title(sprintf('Cluster %s', pclID{trnSubs_psth(csp)}))
end
legend({'Auto-correlogram','Cluster PSTH'})
%% Saving the TRN classification
clInfo = addvars(clInfo, false(size(clInfo,1),1), 'NewVariableNames', 'TRN');
clInfo{gclID(trnSubs),'TRN'} = true;
writeClusterInfo(clInfo,fullfile(dataDir,'cluster_info.tsv'));
%% Plotting the 
% Placing the TRN channels on the probe
probeFig = figure('Name','Probe','Color',[1,1,1]); 
probAx = axes('Parent', probeFig);
plot(chPos(:,1), chPos(:,2),'LineStyle', 'none', 'Marker', 'o',...
    'MarkerEdgeColor', ones(1,3)*0.8)
[~,trnChanSub] = find(chMap' == clInfo{clInfo.TRN==1,'ch'});
set(probAx,'Box','off'); set(get(probAx,'XAxis'),'Color','none');
probAx.YAxis.Color = 'none';
hold(probAx, 'on'); scatter(probAx, chPos(trnChanSub,1),...
    chPos(trnChanSub,2), 'MarkerFaceColor', [0,0.6,0],...
    'MarkerEdgeColor', 'none')
% Placing the electrode name on the probe
text(35:250:(250*3)+35, -20*ones(4,1),...
    arrayfun(@(x) sprintf('Shank %d',x), (1:4)', 'UniformOutput', 0),...
    'HorizontalAlignment', 'center', 'Color', [0.7,0.7,0.7])
probLeg = legend({'E1-channels','TRN cells'' locations'});
probLeg.Box = 'off'; probLeg.Color = 'none'; probLeg.Orientation = 'horizontal';
probLeg.Location = 'north';
probeFig = configureFigureToPDF(probeFig);
probFile = fullfile(figureDir, 'Probemap TRN cluster''s position');
if ~exist([probFile,'.fig'],'file')
    savefig(probeFig, [probFile, '.fig'])
end
if ~exist([probFile,'.pdf'],'file')
    print(probeFig, [probFile, '.pdf'],'-dpdf','-fillpage')
end
if ~exist([probFile,'.emf'],'file')
    print(probeFig, [probFile, '.emf'],'-dmeta')
end

