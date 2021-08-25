%% VPM Control
ClustersControlOnset1stVPMLong=nan(size(ClustersControlOnsetVPM,1),size(ClustersControlOnsetVPM,2));
id=find(ClustersControlOnsetVPM>=-0.2 & ClustersControlOnsetVPM<=0.05);
for i=1:size(id);
ClustersControlOnset1stVPMLong(id(i))=ClustersControlOnsetVPM(id(i));
end
zRespControl=zscore(hist(ClustersControlOnset1stVPMLong,500));
sigRespControl=zRespControl>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigResp=nan(250,size(ClustersControlOnset1stVPM,2));
for i = 1:size(sigRespControl,2)
    TsigResp=Tscale(sigRespControl(:,i)==1);
    TWsigResp(1:size(TsigResp,2),i)=TsigResp;
end
TWsigResp=TWsigResp+0.0025

Ctl=nan(size(TWsigResp,1),size(TWsigResp,2));
for i=1:size(TWsigResp,2)
    BCtl=TWsigResp>0.005&TWsigResp<0.04;
    cCtl=TWsigResp(BCtl(:,i),i);
    Ctl(1:size(cCtl,1),i)=cCtl;
end
    
%% VPM AI
ClustersL200Onset1stVPMLong=nan(size(ClustersL200OnsetVPM,1),size(ClustersL200OnsetVPM,2));
idAI=find(ClustersL200OnsetVPM>=-0.2 & ClustersL200OnsetVPM<=0.05);
for i=1:size(idAI);
ClustersL200Onset1stVPMLong(idAI(i))=ClustersL200OnsetVPM(idAI(i));
end
zRespAI=zscore(hist(ClustersL200Onset1stVPMLong,500));
sigRespAI=zRespAI>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespAI=nan(250,size(ClustersL200Onset1stVPM,2));
for i = 1:size(sigRespAI,2)
    TsigRespAI=Tscale(sigRespAI(:,i)==1);
    TWsigRespAI(1:size(TsigRespAI,2),i)=TsigRespAI;
end
TWsigRespAI=TWsigRespAI+0.0025

C=nan(size(TWsigRespAI,1),size(TWsigRespAI,2));
for i=1:size(TWsigRespAI,2)
    B=TWsigRespAI>0.005&TWsigRespAI<0.04;
    c=TWsigRespAI(B(:,i),i);
    C(1:size(c,1),i)=c;
end
%% S1 Control
ClustersControlOnset1stS1Long=nan(size(ClustersControlOnsetS1,1),size(ClustersControlOnsetS1,2));
idS1=find(ClustersControlOnsetS1>=-0.2 & ClustersControlOnsetS1<=0.05);
for i=1:size(idS1);
ClustersControlOnset1stS1Long(idS1(i))=ClustersControlOnsetS1(idS1(i));
end
zRespControlS1=zscore(hist(ClustersControlOnset1stS1Long,500));
sigRespControlS1=zRespControlS1>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespS1=nan(250,size(ClustersControlOnset1stS1,2));

for i = 1:size(sigRespControlS1,2)
    TsigRespS1=Tscale(sigRespControlS1(:,i)==1);
    TWsigRespS1(1:size(TsigRespS1,2),i)=TsigRespS1;
end
TWsigRespS1=TWsigRespS1+0.0025

CtlS1=nan(size(TWsigRespS1,1),size(TWsigRespS1,2));
for i=1:size(TWsigRespS1,2)
    BCtlS1=TWsigRespS1>0.005&TWsigRespS1<0.04;
    cCtlS1=TWsigRespS1(BCtlS1(:,i),i);
    CtlS1(1:size(cCtlS1,1),i)=cCtlS1;
end
    
%% S1 Laser
ClustersL200Onset1stS1Long=nan(size(ClustersL200OnsetS1,1),size(ClustersL200OnsetS1,2));
idAIS1=find(ClustersL200OnsetS1>=-0.2 & ClustersL200OnsetS1<=0.05);
for i=1:size(idAIS1);
ClustersL200Onset1stS1Long(idAIS1(i))=ClustersL200OnsetS1(idAIS1(i));
end
zRespAIS1=zscore(hist(ClustersL200Onset1stS1Long,500));
sigRespAIS1=zRespAIS1>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespAIS1=nan(250,size(ClustersL200Onset1stS1,2));
for i = 1:size(sigRespAIS1,2)
    TsigRespAIS1=Tscale(sigRespAIS1(:,i)==1);
    TWsigRespAIS1(1:size(TsigRespAIS1,2),i)=TsigRespAIS1;
end
TWsigRespAIS1=TWsigRespAIS1+0.0025

CS1=nan(size(TWsigRespAIS1,1),size(TWsigRespAIS1,2));
for i=1:size(TWsigRespAIS1,2)
    BS1=TWsigRespAIS1>0.005&TWsigRespAIS1<0.04;
    cS1=TWsigRespAIS1(BS1(:,i),i);
    CS1(1:size(cS1,1),i)=cS1;
end
%% TRN Control
ClustersControlOnset1stTRNLong=nan(size(ClustersControlOnsetTRN,1),size(ClustersControlOnsetTRN,2));
idTRN=find(ClustersControlOnsetTRN>=-0.2 & ClustersControlOnsetTRN<=0.05);
for i=1:size(idTRN);
ClustersControlOnset1stTRNLong(idTRN(i))=ClustersControlOnsetTRN(idTRN(i));
end
zRespControlTRN=zscore(hist(ClustersControlOnset1stTRNLong,500));
sigRespControlTRN=zRespControlTRN>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespTRN=nan(250,size(ClustersControlOnset1stTRN,2));
for i = 1:size(sigRespControlTRN,2)
    TsigRespTRN=Tscale(sigRespControlTRN(:,i)==1);
    TWsigRespTRN(1:size(TsigRespTRN,2),i)=TsigRespTRN;
end
TWsigRespTRN=TWsigRespTRN+0.0025

CtlTRN=nan(size(TWsigRespTRN,1),size(TWsigRespTRN,2));
for i=1:size(TWsigRespTRN,2)
    BCtlTRN=TWsigRespTRN>0.005&TWsigRespTRN<0.04;
    cCtlTRN=TWsigRespTRN(BCtlTRN(:,i),i);
    CtlTRN(1:size(cCtlTRN,1),i)=cCtlTRN;
end
    
%% TRN Laser
ClustersL200Onset1stTRNLong=nan(size(ClustersL200OnsetTRN,1),size(ClustersL200OnsetTRN,2));
idAITRN=find(ClustersL200OnsetTRN>=-0.2 & ClustersL200OnsetTRN<=0.05);
for i=1:size(idAITRN);
ClustersL200Onset1stTRNLong(idAITRN(i))=ClustersL200OnsetTRN(idAITRN(i));
end
zRespAITRN=zscore(hist(ClustersL200Onset1stTRNLong,500));
sigRespAITRN=zRespAITRN>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespAITRN=nan(250,size(ClustersL200Onset1stTRN,2));
for i = 1:size(sigRespAITRN,2)
    TsigRespAITRN=Tscale(sigRespAITRN(:,i)==1);
    TWsigRespAITRN(1:size(TsigRespAITRN,2),i)=TsigRespAITRN;
end
TWsigRespAITRN=TWsigRespAITRN+0.0025;

CTRN=nan(size(TWsigRespAITRN,1),size(TWsigRespAITRN,2));
for i=1:size(TWsigRespAITRN,2)
    BTRN=TWsigRespAITRN>0.005&TWsigRespAITRN<0.04;
    cTRN=TWsigRespAITRN(BTRN(:,i),i);
    CTRN(1:size(cTRN,1),i)=cTRN;
end

