fullName = @(x) string(fullfile(x.folder, x.name));
fnOpts = {'UniformOutput', false};
%%
expDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch2_ephys\MC\WT13\211119_C";
behDirs = dir(fullfile(expDir, "Beh*")); behDir = fullName(behDirs);
%%
ephDirs = dir(fullfile(expDir, "ephys_*")); ephDir = fullName(ephDirs);
condFiles = dir(fullfile(ephDir, "*analysis.mat"));
load(fullName(condFiles), "Conditions", "fs")
%%
ccSub = find(arrayfun(@(c) contains(c.name, ["Control Puff", "PTX","Death"], ...
    "IgnoreCase", true), Conditions));

pairedStim = arrayfun(@(x) Conditions(1).Triggers(:,1) == ...
    Conditions(x).Triggers(:,1)', ccSub, fnOpts{:});
pairedStim = cellfun(@(x) any(x, 2), pairedStim, fnOpts{:});
pairedStim = cat(2, pairedStim{:});
consCondNames = arrayfun(@(x) string(x.name), Conditions(ccSub));
%%
[behRes, behFigDir] = analyseBehaviour(behDir, 'PairedFlags', pairedStim, ...
    'ConditionsNames', cellstr(consCondNames));
biFigPttrn = "BehIndex%s";
biFigPttrn = sprintf(biFigPttrn, sprintf(" %s (%%.3f)", consCondNames));
[pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
behRes = arrayfun(@(bs, ba) setfield(bs,'BehIndex', ba), behRes, pAreas);
set(behAreaFig, 'UserData', behRes)
biFN = sprintf(biFigPttrn, pAreas);
saveFigure(behAreaFig, fullfile(behFigDir, biFN), true);