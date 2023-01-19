%%struct evoked firing rates MOCK
structString
MIevokE=(evFr(:,2) - evFr(:,1))./(evFr(:,1) + evFr(:,2))
MIevokESP=(spFr(:,2) - spFr(:,1))./(spFr(:,1) + spFr(:,2))
posMod=MIevokE>=0.001
negMod=MIevokE<=-0.001
VPMMockEvokr=evFr(wruIdx&signMod==0,:)
VPMMockEvokPot=evFr(wruIdx&posMod&signMod,:)
VPMMockEvokDep=evFr(wruIdx&negMod&signMod,:)
vpmMockEVStrtR=struct('NonMod',VPMMockEvokr,'rPot',VPMMockEvokPot,'rDep',VPMMockEvokDep)
vpmMockEVStrtR2=vpmMockEVStrtR
%%struct spontaneous firing rates MOCK
VPMMockSpontr=pfr(wruIdx&signMod==0,:)
VPMMockSpontDep=pfr(wruIdx&negMod&signMod,:)
VPMMockSpontPot=pfr(wruIdx&posMod&signMod,:)
vpmMockSPstrR2=struct('NonMod',VPMMockSpontr,'rPot',VPMMockSpontPot,'rDep',VPMMockSpontDep)
vpmMockSPstrR=vpmMockSPstrR2
%% struct evoked firing rate VPM-L6

MIevokE=(evFr(:,2) - evFr(:,1))./(evFr(:,1) + evFr(:,2))
MIevokSP=(spFr(:,2) - spFr(:,1))./(spFr(:,1) + spFr(:,2))
posMod=MIevokSP>=0.001
negMod=MIevokSP<=-0.001
VPMEvokr=evFr(wruIdx&signMod==0,:)
VPMEvokPot=evFr(wruIdx&posMod&signMod,:)
VPMEvokDep=evFr(wruIdx&negMod&signMod,:)
vpmEVStrtR2=struct('NonMod',VPMEvokr,'rPot',VPMEvokPot,'rDep',VPMEvokDep)
%%struct spontaneous firing rates MOCK
VPMSpontr=spFr(wruIdx&signMod==0,[1 2])
VPMSpontDep=spFr(wruIdx&negMod&signMod,[1 2])
VPMSpontPot=spFr(wruIdx&posMod&signMod,[1 2])
vpmSPstrR=struct('NonMod',VPMSpontr,'rPot',VPMSpontPot,'rDep',VPMSpontDep)
vpmSPstrR=vpmSPstrR2