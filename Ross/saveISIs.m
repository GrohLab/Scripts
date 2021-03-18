close all
isiDir = fullfile(figureDir,'ISI\');
if ~mkdir(isiDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end

plotISI(spkSubs2,fs,gclID);
for a = 1:length(gclID(2:end))
    H = figure(a);
    savefig(H,fullfile(isiDir, ['Unit ', num2str(gclID{a}), '.fig']), 'compact');
end