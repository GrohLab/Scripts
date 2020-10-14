function Venn(clInfo, a)
% a = laser intensity
pwr = num2str(a);
MRind = ['Mech_Control_',pwr, 'mW_MR'];
MR = sum(clInfo.(MRind));
LRind = ['Laser_Control_', pwr, 'mW_LR'];
LR = sum(clInfo.(LRind));
LR_and_MR = sum(clInfo.(MRind) & clInfo.(LRind));
deltaMRind = ['Mech_Control_' pwr, 'mW_vs_Mech_Laser_', pwr, 'mW_Evoked_Response'];
deltaMR = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind));
modMR = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind) == false);
vennX([LR - LR_and_MR, LR_and_MR, MR - LR_and_MR, modMR, 0, 0, deltaMR], 0.01);


LEvCounts = ['Mech_Laser_' pwr, 'mW_Counts_Evoked'];
CEvCounts = ['Mech_Control_' pwr, 'mW_Counts_Evoked'];

modInc = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind) == false & clInfo.(LEvCounts) > clInfo.(CEvCounts));
modDec = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind) == false & clInfo.(LEvCounts) < clInfo.(CEvCounts));
drvInc = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind) & clInfo.(LEvCounts) > clInfo.(CEvCounts));
drvDec = sum(clInfo.(MRind) & clInfo.(deltaMRind) & clInfo.(LRind) & clInfo.(LEvCounts) < clInfo.(CEvCounts));
vennX([modInc, drvInc, LR, drvDec, modDec, 0, 0], 0.005);
end
