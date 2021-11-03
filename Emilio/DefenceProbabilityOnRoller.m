% No previous movement before nor after the stimulus
excludeFlag = rms(vStack(spontFlag,:)) > 0.85 | rms(vStack) > 0.9;
% Previous movement before 
% excludeFlag = rms(vStack(spontFlag,:)) < 0.85;
Na = sum(delayFlags & ~excludeFlag(:));
rngRollSpeed = cell(Nccond,1);
responseWindow = [0, 0.4];
responseFlags = stTx >= responseWindow(1) & stTx <= responseWindow(2);
for ccond = 1:Nccond
    rngRollSpeed{ccond} =...
        max(abs(vStack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end
% Maximum number of trials
maxNt = max(cellfun(@numel, rngRollSpeed));
% Placing the trial's maximum speed in a matrix
speedsMat = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    rngRollSpeed, fnOpts{:}); speedsMat = cat(2, speedsMat{:});
speedsMat = speedsMat * en2cm;
% Cutting on different thresholds for all trials
thetaSpeed = 0.2:0.2:6;
moveFlag = speedsMat > reshape(thetaSpeed,1,1,[]);
probMove = sum(moveFlag)./Na; probMove = squeeze(probMove);
probFig = figure; plot(thetaSpeed, probMove'); ylim([0,1])
set(gca, "Box", "off", "Color", "none")
consCondNames = arrayfun(@(x) x.name, Conditions(consideredConditions),fnOpts{:});
lgnd = legend(consCondNames); set(lgnd, "Box", "off", "Location", "best")
ylabel("Trial proportion / Movement probability")
xlabel("Speed threshold \theta [cm/s]")
title("Trial proportion with elicited movement")
%% Save?
sveFigAns = questdlg("Save figure?", "Save", "Yes", "No", "Yes");
if strcmpi(sveFigAns, "Yes")
    saveFigure(probFig, fullfile(figureDir, 'Movement probability'), 1, 1)
end
%% Bar plotting
condSubs = [1:7]; thS = 20;
% Table with the successful trials in the first column.
tbl = [probMove(condSubs, thS).*Na(condSubs)',...
    abs(probMove(condSubs, thS).*Na(condSubs)' - Na(condSubs)')];
% Chi squared test
sumCol = sum(tbl)'; sumFil = sum(tbl,2)';
multSum = (sumCol * sumFil)'; esperado = multSum / sum(tbl(:));
resChi = ((tbl - esperado).^2)./esperado; resChi = sum(resChi(:));
pChi = chi2pdf(resChi,1); pFih = Inf;
if all(size(tbl) == [2, 2])
    [~, pFih] = fishertest(tbl);
end
% Plotting the proportions and test results
barFig = figure; bar(probMove(:,thS), "EdgeColor", "none");
yHeight = max(probMove(condSubs, thS))*1.05;
hold on; plot(condSubs, yHeight * ones(size(condSubs,2),1),...
    'LineWidth', 1, "Color", 'k')
signStr = 'n.s.';
if pChi < 5e-3 || pFih < 5e-3
    signStr = '***';
elseif pChi < 1e-2 || pFih < 1e-2
    signStr = '**';
elseif pChi < 5e-2 || pFih < 5e-2
    signStr = '*';
end
text(mean(condSubs), yHeight, signStr, "HorizontalAlignment", "center",...
    "VerticalAlignment", "bottom");
xticklabels(consCondNames); ylim([0,(yHeight >= 1)*1.1 + (yHeight < 1)*1]);
set(get(barFig,"Children"), "Box", "off", "Color", "none")
title(sprintf('Movement probability \\theta:%.1f cm/s', thetaSpeed(thS)))
ylabel("Trial proportion / Movement probability")
%% Save?
sveFigAns = questdlg("Save figure?", "Save", "Yes", "No", "Yes");
if strcmpi(sveFigAns, "Yes")
    saveFigure(barFig, fullfile(figureDir,...
        sprintf('Movement probability bar TH%.1f cmps', thetaSpeed(thS))),...
        1, 1)
end

%%
xticklabels({'P+L','P','L'})