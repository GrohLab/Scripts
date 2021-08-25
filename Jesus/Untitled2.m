ProbabilitybinAI=nan(50,size(ClustersL200Onset1stVPM,2));
for i= 1:size(ClustersL200Onset1stVPM,2);
    latAI=histcounts(ClustersL200Onset1stVPM(ClustersL200Onset1stVPM(:,i)>0.003),50,'Normalization','probability');
    ProbabilitybinAI(:,i)=latAI;
end
%%    
minlatencyAI=nan(2000,size(ClustersL200Onset1stVPM,2));
for i= 1:size(ClustersL200Onset1stVPM,2);
    latenciesAI=ClustersL200Onset1stVPM(ClustersL200Onset1stVPM(:,i)>0.005 & ClustersL200Onset1stVPM(:,i)<0.02);
    minlatencyAI(1:size(latenciesAI,1),i) = latenciesAI;
    WminlatencyAI=min(minlatencyAI,[],1)
end
%%
minlatencyS1=nan(size(ClustersControlOnset1stS1,1),size(ClustersControlOnset1stS1,2));
for i= 1:size(ClustersControlOnset1stS1,2);
    siglatColS1=ClustersControlOnset1stS1>0.004;
     BS1=ClustersControlOnset1stS1(siglatColS1(:,i),i)
    minlatencyS1(1:size(BS1,1),i) = BS1;
end
%%
minlatencyMOCK=nan(size(ClustersControlOnset1stMOCK,1),size(ClustersControlOnset1stMOCK,2));
for i= 1:size(ClustersControlOnset1stMOCK,2);
    siglatColMOCK=ClustersControlOnset1stMOCK>0.004;
     BS1MOCK=ClustersControlOnset1stMOCK(siglatColMOCK(:,i),i)
    minlatencyMOCK(1:size(BS1MOCK,1),i) = BS1MOCK;
end
%%
minlatencyAIMOCK=nan(size(ClustersL200Onset1stMOCK,1),size(ClustersL200Onset1stMOCK,2));
for i= 1:size(ClustersL200Onset1stMOCK,2);
    siglatColAIMOCK=ClustersL200Onset1stMOCK>0.004;
     BS1AIMOCK=ClustersL200Onset1stMOCK(siglatColAIMOCK(:,i),i)
    minlatencyAIMOCK(1:size(BS1AIMOCK,1),i) = BS1AIMOCK;
end
%%
minlatencyAI=nan(20000,276);
for i= 1:size(ClustersL200Onset1stVPM,2);
    latCol=ClustersL200Onset1stVPM>0.005 & ClustersL200Onset1stVPM<0.015;
    B=ClustersL200Onset1stVPM(latCol(:,i));
    minlatencyAI(1:size(B,1),i) = B;
end
%%
minlatencyAI=nan(size(ClustersL200Onset1stVPM,1),size(ClustersL200Onset1stVPM,2));
idAI=find(ClustersL200Onset1stVPM>=0.005 & ClustersL200Onset1stVPM<=0.015);
for i=1:size(idAI);
    minlatencyAI(idAI(i))=ClustersL200Onset1stVPM(idAI(i));
end
%%
minlatency=nan(size(ClustersControlOnset1stVPM,1),size(ClustersControlOnset1stVPM,2));
idControl=find(ClustersControlOnset1stVPM>=0.005 & ClustersControlOnset1stVPM<=0.015);
for i=1:size(idControl);
    minlatency(idControl(i))=ClustersControlOnset1stVPM(idControl(i));
end
%%
minlatencyALL=nan(size(Clusters1stVPM,1),size(Clusters1stVPM,2));
idALL=find(Clusters1stVPM>=0.005 & Clusters1stVPM<=0.015);
for i=1:size(idALL);
    minlatencyALL(idALL(i))=Clusters1stVPM(idALL(i));
end
%%
ClustersControlOnset1stVPMLong=nan(size(ClustersControlOnsetVPM,1),size(ClustersControlOnsetVPM,2));
id=find(ClustersControlOnsetVPM>=-0.2 & ClustersControlOnsetVPM<=0.05);
for i=1:size(id);
ClustersControlOnset1stVPMLong(id(i))=ClustersControlOnsetVPM(id(i));
end
zRespControl=zscore(hist(ClustersControlOnset1stVPMLong,500));
sigRespControl=zRespControl>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigResp=nan(250,276);
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
    
%%
ClustersL200Onset1stVPMLong=nan(size(ClustersL200OnsetVPM,1),size(ClustersL200OnsetVPM,2));
idAI=find(ClustersL200OnsetVPM>=-0.2 & ClustersL200OnsetVPM<=0.05);
for i=1:size(idAI);
ClustersL200Onset1stVPMLong(idAI(i))=ClustersL200OnsetVPM(idAI(i));
end
zRespAI=zscore(hist(ClustersL200Onset1stVPMLong,500));
sigRespAI=zRespAI>=1.95;
Tscale=[-0.200:0.0005:0.05-0.0005];

TWsigRespAI=nan(250,276);
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

    

%%
minlatency8=[];
for i= 1:size(ClustersControlOnset1stS1,2);
    latVPM8=min(ClustersControlOnset1stVPM8(ClustersControlOnset1stVPM8(:,i)>0.004));
    minlatency8(i,1)=latVPM8;
end
