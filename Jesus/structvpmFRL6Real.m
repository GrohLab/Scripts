%structure evoked firing rates VPM
structString
posMod=MIevok>=0.001
negMod=MIevok<=-0.001
VPML6Evokr=evFr(wruIdx&signMod==0,:)
VPML6EvokPot=evFr(wruIdx&posMod&signMod,:)
VPML6EvokDep=evFr(wruIdx&negMod&signMod,:)
vpmL6EVStrtR2=struct('NonMOd',VPML6Evokr,'rPot',VPML6EvokPot,'rDep',VPML6EvokDep)
%struct spontaneous firing rates VPM
VPML6Spontr=pfr(wruIdx&signMod==0,:)
VPML6SpontDep=pfr(wruIdx&negMod&signMod,:)
VPML6SpontPot=pfr(wruIdx&posMod&signMod,:)
vpmL6SPstrR2=struct('NonMod',VPML6Spontr,'rPot',VPML6SpontPot,'rDep',VPML6SpontDep)
structString