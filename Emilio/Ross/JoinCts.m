CFA.SpontCts = CFAisi(1).Vals(1).cts;
for a = 2:length(CFAisi)
CFA.SpontCts = CFA.SpontCts + CFAisi(a).Vals(1).cts;
end
CFA.EvokedCts = CFAisi(1).Vals(2).cts;
for a = [4, 7]
CFA.EvokedCts = CFA.EvokedCts + CFAisi(a).Vals(2).cts;
end
CFA.bns = CFAisi(1).Vals(2).bns;
save(fullfile('Z:\Ross\Experiments\VPL\VPL_CFA','CFA_Spont&Evoked.mat'), 'CFA', '-v7.3');