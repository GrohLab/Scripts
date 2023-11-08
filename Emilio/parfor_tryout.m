Npr = 1000;
for cl = 1:find(wruIdx)
    auxSubs = arrayfun(@(x) randperm(NnzvPcl(cl)), 1:Npr, 'UniformOutput', false);
    auxSpks = cellfun(@(s) round(cumsum(ISIVals{cl}(s))*fs), auxSubs, 'UniformOutput', false);
    rndStack = getStacks(false, Conditions(chCond).Triggers, onOffStr,...
        timeLapse, fs, fs, cat(2, spkSubs{cl}, auxSpks));
    rndStack(2,:,:) = [];
    rst2 = arrayfun(@(x) getRasterFromStack(rndStack, ~delayFlags(:,x), ...
        [false; true(Npr-2,1)], timeLapse, fs, true, true), ...
        1:size(delayFlags,2), 'UniformOutput', false);
    [C, bE] = arrayfun(@(c) histcounts([rst2{1}{c,:}], 'BinWidth', ...
        binSz, 'BinLimits', timeLapse), 1:size(rst2{1},1), ...
        'UniformOutput', false);
    bE = bE{1}; Ctot = cat(1, C{:}); bC = mean([bE(1:end-1);bE(2:end)]);
    P = arrayfun(@(x) fitdist(Ctot(2:end,x), 'Poisson'), 1:size(Ctot,2));
    lambdas = arrayfun(@(x) x.lambda, P);
    PI = arrayfun(@(x) x.paramci, P, 'UniformOutput', false); PI = cat(2, PI{:});
    figure; plot(bC, Ctot')
end