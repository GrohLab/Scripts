function Venn(clInfo, a)
% a = laser intensity
pwr = num2str(a);

LRind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) & clInfo.(['Mech_Control_', pwr, 'mW_MR']) == false...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']) == false);
LR = length(LRind);

MRind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) == false & clInfo.(['Mech_Control_', pwr, 'mW_MR'])...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']) == false);
MR = length(MRind);

LRMRind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) & clInfo.(['Mech_Control_', pwr, 'mW_MR'])...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']) == false);
LRMR = length(LRMRind);

Drivenind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) & clInfo.(['Mech_Control_', pwr, 'mW_MR'])...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']));
Driven = length(Drivenind);

ModulatedMRind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) == false & clInfo.(['Mech_Control_', pwr, 'mW_MR'])...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']));
ModulatedMR = length(ModulatedMRind);

ModulatedLRind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) & clInfo.(['Mech_Control_', pwr, 'mW_MR']) == false...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']));
ModulatedLR = length(ModulatedLRind);

Awokenind = find(clInfo.(['Laser_Control_', pwr, 'mW_LR']) == false & clInfo.(['Mech_Control_', pwr, 'mW_MR']) == false...
& clInfo.(['Mech_Control_', pwr, 'mW_vs_' 'Mech_Laser_', pwr, 'mW_Evoked_Response']));
Awoken = length(Awokenind);

vennX([LR, LRMR, MR, ModulatedMR , Awoken, ModulatedLR, Driven], 0.01);

LEvCounts = ['Mech_Laser_' pwr, 'mW_Counts_Evoked'];
CEvCounts = ['Mech_Control_' pwr, 'mW_Counts_Evoked'];

modInc = sum(clInfo.(LEvCounts)(ModulatedMRind) > clInfo.(CEvCounts)(ModulatedMRind));
modDec = sum(clInfo.(LEvCounts)(ModulatedMRind) < clInfo.(CEvCounts)(ModulatedMRind));
drvInc = sum(clInfo.(LEvCounts)(Drivenind) > clInfo.(CEvCounts)(Drivenind));
drvDec = sum(clInfo.(LEvCounts)(Drivenind) < clInfo.(CEvCounts)(Drivenind));
Laser_Rest = sum(clInfo.(['Laser_Control_', pwr, 'mW_LR'])) - drvInc - drvDec;
vennX([modInc, drvInc, Laser_Rest, drvDec, modDec, 0, 0], 0.01);
end
