fnOpts = {"UniformOutput", false};
cpts = getWaveformCriticalPoints(vStack, fsRoll);
cpts = cellfun(@(x) x - 1, cpts, 'UniformOutput', false);
camp = arrayfun(@(x) interp1(stTx, vStack(:,x), cpts{x,1}),...
    (1:size(vStack,2))', "UniformOutput", false);
[~, mxSub] = cellfun(@max, camp, 'UniformOutput', false);
dlTms = arrayfun(@(x) cpts{x,1}(mxSub{x}), (1:size(cpts,1))',"UniformOutput",false);
ttms = sort(cat(1,Conditions(3:5).Triggers))./fs;
emptyCells = cellfun(@isempty, dlTms);
pts = [cat(1,dlTms{~emptyCells}), ttms(~emptyCells)];
