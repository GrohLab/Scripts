%%
Nat = size(atTimes{1},1); Nit = size(itTimes{1}(:,1),1);
err = Inf; r2 = 0;
errTh = -3;
atDelta = diff(atTimes{1}(:,1)); itDelta = diff(itTimes{1}(:,1));
ddm = distmatrix(itDelta, atDelta);
%figure;
cit = 1;
% Arduino times always come before the intan.
dm = atTimes{1} - itTimes{1}(:,1)'; preArd = dm < 0;
[~, mnSubs] = sort(dm(:), "ascend", "ComparisonMethod", "abs");
mnSubs(~preArd(mnSubs)) = []; [aSubs, iSubs] = ind2sub(size(dm),  mnSubs);
mxSub = min(Nit, Nat);
cp = 1;
figure; stem(itTimes{1}(:,1), ones(size(itTimes{1},1),1));
hold on; stem(atTimes{1}, -ones(size(atTimes{1},1),1));
while err > errTh || r2 < 0.99
    while length(unique(iSubs(cp:mxSub-1+cp))) < mxSub ||...
            length(unique(aSubs(cp:mxSub-1+cp))) < mxSub
        cp = cp + 1;
    end
    xSubs = sort(iSubs(cp:mxSub-1+cp));
    x = itTimes{1}(xSubs,1);
    ySubs = sort(aSubs(cp:mxSub-1+cp));
    y = atTimes{1}(ySubs);
    [mdl, yhat, r2] = fit_poly(x, y, 1);
    err = log(norm(y - yhat, 1));
    cp = cp + 1;
end
%%
figure; stem(itTimes{cc}(:,1), ones(size(itTimes{cc},1),1));
hold on; stem(atTimes{ca}, -ones(size(atTimes{ca},1),1));
%%
figure; imagesc(log(abs(dm{cc}+1)))
%%
[bC, bE] = prepareLogBinEdges(abs(dm{cc}(preArd)), 32);
bH = histcounts(log10(abs(dm{cc}(preArd))), bE);
bH = bH./(diff(10.^bE) * sum(bH));
[~, delSubs] = sort(bH, "descend");

figure; semilogx(10.^bC, bH)

% Arduino times that are before the intan
posPairs = arrayfun(@(x) find(preArd(:,x), 1, "last"), 1:Nit);

%%
while err > errTh
    % DEFINE: randSubs
    if Nat > Nit
        % Choose random sets of points fitting with intan cardinality
        
        x = itTimes{1}(:,1);
        y = atTimes{1}(randSubs);
    elseif Nat < Nit
        % Choose random sets of points fitting with arduino cardinality
        x = itTimes{1}(setdiff(1:Nit, randperm(Nit, Nit - Nat)),1);
        y = atTimes{1};
    else
        % Compute y = mx + b
        x = itTimes{1}(:,1); y = atTimes{1};
    end
    [mdl, yhat, r2] = fit_poly(x, y, 1);
    hold on; scatter(x, y, ".", "MarkerFaceAlpha", 0.3)
    plot([min(x);max(x)], [min(x), 1; max(x), 1]*mdl,...
        "DisplayName", num2str(r2))
    errVec = abs(y - yhat); err = log(norm(y - yhat, 1));
end

%%
[Nta, Nti] = cellfun(@size, dm); cc = 1; ca = 3 - cc;
preArd = dm{cc} < 0; [~, mnSubs] = sort(dm{cc}(:), "ascend",...
    "ComparisonMethod", "abs");
mnSubs(~preArd(mnSubs)) = [];
[aSubs, iSubs] = ind2sub(size(dm{cc}), mnSubs);
mxSub = min(Nti(cc),Nta(cc)); naiveSubs = 1:mxSub;
pairFlags = abs(iSubs - aSubs) < abs(Nta(cc) - Nti(cc)) + 1;
dstPrs = dm{cc}(mnSubs(pairFlags)); medDst = median(dstPrs);
unqePrs = unique([...
    iSubs(pairFlags), aSubs(pairFlags), dstPrs], "row");
allFlag = unqePrs(:,1:2) == reshape(naiveSubs, 1,1,[]);
IoA = all(any(allFlag,3))'; snglPairs = zeros(mxSub, 2);
if all(IoA)
    [~, IoASub] = min([Nti(cc),Nta(cc)]); IoA = false(2,1);
    IoA(IoASub) = true;
end
auxCnt = 1;
% Logical reduction of possibilities
while auxCnt <= mxSub && ~isempty(unqePrs)
    allFlag = unqePrs(:,1:2) == reshape(naiveSubs, 1,1,[]);
    appCount = squeeze(sum(allFlag, 1));
    appCount(appCount == 0) = NaN;
    % Lower cardinality trigger times with only one pairing possibility
    mstUseFlag = any(appCount == IoA);
    if all(~mstUseFlag)
        unlikelyFlag = isnan(appCount(IoA,:));
        unpairFlag = any(unqePrs(:,IoA) == naiveSubs(unlikelyFlag),2);
        if all(~unpairFlag)
            unlikelyFlag = unqePrs(:,3) < medDst*1.1;
            if all(~unlikelyFlag)
                fprintf(1, "Taking the pair(s) that is(are) closest to")
                fprintf(1, " the linear fit\n")
                xSubs = snglPairs(snglPairs(:,1) ~= 0,1);
                ySubs = snglPairs(snglPairs(:,2) ~= 0,2);
                x = itTimes{cc}(xSubs,1);
                y = atTimes{ca}(ySubs);
                [n, d] = getHesseLineForm(fit_poly(x, y, 1));
                M = [itTimes{cc}(unqePrs(:,1),1), atTimes{ca}(unqePrs(:,2))];
                yErr = log(abs(M*n - d));
                bigErrFlag = yErr > 0;
                unqePrs(bigErrFlag,:) = [];
            else
                fprintf(1, "Removing pairs with a 'big'")
                fprintf(1, " distance. (%.2f)\n", medDst*1.1)
                unqePrs(unlikelyFlag,:) = [];
            end
        else
            unqePrs(unpairFlag,:) = [];
        end
        continue
    end
    % Pairs that *must* be in the pair selection
    unPrFlag = any(unqePrs(:,IoA) == naiveSubs(mstUseFlag),2);
    Nup = sum(unPrFlag);
    snglPairs(auxCnt:auxCnt-1+Nup,:) = unqePrs(unPrFlag,1:2);
    auxCnt = auxCnt+Nup;
    unavFlags = arrayfun(@(x) any(unqePrs(:,x) == snglPairs(:,x)',2), 1:2,...
        "UniformOutput", false); unavFlags = [unavFlags{:}];
    unqePrs(any(unavFlags,2),:) = [];
end
%%
foldersInBatch = dir("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch2\*\*\*");
fileFlags = [foldersInBatch.isdir]';
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), foldersInBatch);
und_Flag = arrayfun(@(x) contains(x.name, "_"), foldersInBatch);
foldersInBatch(pointFlag | ~fileFlags | ~und_Flag) = [];
%% 
for cf = foldersInBatch'
    dataDir = fullfile(cf.folder, cf.name);
    fprintf(1, "Processing folder %s\n", dataDir)
    readAndCorrectArdTrigs(dataDir);
end
%%
pmFlags = arrayfun(@(x) ismember(Conditions(pFlag).Triggers(:,1),...
    Conditions(x).Triggers(:,1)), mixSubs, fnOpts{:}); pmFlags = [pmFlags{:}];


lmFlags = arrayfun(@(x) ismember(Conditions(lFlag).Triggers(:,1),...
    Conditions(x).Triggers(:,1)), mixSubs, fnOpts{:}); lmFlags = [lmFlags{:}];
for ctf = 1:Nat
    atName = fullfile(atFiles(ctf).folder, atFiles(ctf).name);
    % Loading only arduino triggers (atTimes and atNames)
    auxVars = load(atName, "at*");
    if isempty(auxVars)
        fprintf(1, "Warning! Unable to load necessary variables")
        fprintf(1, " from the condition file:\n%s\n", atFiles(ctf).name)
        return
    end
    fieldnames(auxVars);
end
%%

figure; subplot(2,3,1); imagesc(psthTx, [], wPSTH(find(anyH)+1,:)./max(wPSTH(find(anyH)+1,:),[],2));
title("Whisker deflection"); yticks(1:sum(anyH));yticklabels(gclID(anyH))
subplot(2,3,2); imagesc(psthTx, [], PSTH(anyH,:,1)./max(PSTH(anyH,:,1),[],2));
title("L5-BC"); set(get(gca, "YAxis"), "Visible", "off")
subplot(2,3,3); imagesc(psthTx, [], PSTH(anyH,:,2)./max(PSTH(anyH,:,2),[],2));
title("L5-MC"); set(get(gca, "YAxis"), "Visible", "off")
% Non-normalized
subplot(2,3,4); imagesc(psthTx, [], wPSTH(find(anyH)+1,:));
title("Whisker deflection"); yticks(1:sum(anyH));yticklabels(gclID(anyH))
subplot(2,3,5); imagesc(psthTx, [], PSTH(anyH,:,1));
title("L5-BC"); set(get(gca, "YAxis"), "Visible", "off")
subplot(2,3,6); imagesc(psthTx, [], PSTH(anyH,:,2));
title("L5-MC"); set(get(gca, "YAxis"), "Visible", "off")

%% 
for cif = 1:size(imFiles,1)
    ifName = fullfile(imFiles(cif).folder, imFiles(cif).name);
    [~, bName] = fileparts(ifName);
    im = imread(ifName, "jpg");
    imInv = 255 - im;
    imwrite(imInv, fullfile(imFiles(cif).folder, string(bName) + "_Inv.jpg"))
end
%% Pooled anaesthetised experiments
dataDir = 'Z:\Emilio\SuperiorColliculusExperiments\Anaesthetised';
load(fullfile(dataDir, 'Opto-responses.mat'))
fnOpts = {"UniformOutput", false}; H = logical(H);
gclID = clInfo(clInfo.ActiveUnit == 1, :).Properties.RowNames;
respWin = [3, 35]*1e-3;
condNames = arrayfun(@(x) string(x.name), relativeSpkTmsStruct);
concatCellRows = @(cmat) arrayfun(@(x) cat(2, cmat{x,:}),...
    (1:size(cmat,1))', fnOpts{:});
getSpikeCellMat = @(spkStruct) arrayfun(@(x) x.SpikeTimes,...
    spkStruct, fnOpts{:});
sqzSpks = cellfun(concatCellRows, getSpikeCellMat(relativeSpkTmsStruct),...
    fnOpts{:});
respSpks = cellfun(@(x) cellfun(@(y) y(y>=respWin(1)&y<=respWin(2)),...
    x, fnOpts{:}), sqzSpks, fnOpts{:});
stDv = cellfun(@(x) cellfun(@(y) std(y), x), respSpks, fnOpts{:});
stDv = cat(2, stDv{:});
[PSTH, bnEdgs] = cellfun(@(x) cellfun(@(y) histcounts(y,...
    "BinLimits", respWin, "BinWidth", 5e-4), x, fnOpts{:}),...
    respSpks, fnOpts{:}); bnEdgs = bnEdgs{1}{1}; 
bnCnts = (bnEdgs(1:end-1)+bnEdgs(2:end))/2;
PSTH = cellfun(@(x) cat(1, x{:}), PSTH, fnOpts{:});
for ccond = 1:size(condNames, 2)
    fig = figure; axs(1) = subplot(1,2,1); 
    imagesc(bnCnts*1e3, [], PSTH{ccond}(H(:,ccond),:))
    yticks(1:sum(H(:,ccond))); title(condNames(ccond))
    yticklabels(cellfun(@(x) strrep(x, '_','\_'), gclID(H(:,ccond)),...
        fnOpts{:})); xlabel("Time [ms]"); axs(2) = subplot(1,2,2); 
    barh(stDv(H(:,ccond),ccond), "EdgeColor", "none"); xlabel("Time [ms]");
    set(get(axs(2), "YAxis"), "Direction", "reverse", "Visible", "off")
    xticklabels(xticks*1e3); title("\sigma in ms")
    arrayfun(@(x) set(x, "Box", "off", "Color", "none"), axs)
    linkaxes(axs,'y')
    saveFigure(fig, fullfile(dataDir,"PopFigures",condNames(ccond)),1,1);
end