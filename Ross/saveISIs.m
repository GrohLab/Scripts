close all
isiDir = fullfile(figureDir,'ISI\');
if ~mkdir(isiDir)
    fprintf(1,'There was an issue with the figure folder...\n');
end

plotISI(spkSubs,fs,gclID);
for a = 1:length(gclID(2:end))
    H = figure(a);
    savefig(H,fullfile(isiDir, ['Cluster_', gclID{a + 1,1}]), 'compact');
end