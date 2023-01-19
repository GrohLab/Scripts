%% statistics
%fisher exact test and tabular table
PotT=table([90;21],[355;233],'RowNames',{'L6','Mock'},'VariableNames',{'Pot','NoPot'})
XModNoMod=table([101;71],[217;150],'RowNames',{'L6','Mock'},'VariableNames',{'Mod','NoMod'})
fishertest(xPotNoPot)

%%chi2test manual

observed=[PotT.Pot(1) PotT.NoPot(1) PotT.Pot(2) PotT.NoPot(2)]
expected=[]
bigTotal=sum(sum(PotT{:,:}))
    
    %% evoked
%Control
    %l6
    FRcontEVL6=nan(318,3)
    FRcontEVL6(1:size(vpmL6Evok.evokAll(:,1),1),1) = vpmL6Evok.evokAll(:,1)
    FRcontEVL6(1:size(vpmL6Evok.rPot(:,1),1),2) = vpmL6Evok.rPot(:,1)
    FRcontEVL6(1:size(vpmL6Evok.rDep(:,1),1),3) = vpmL6Evok.rDep(:,1)
    %mock
    FRcontEVMock=nan(318,3)
    FRcontEVMock(1:size(vpmMockEVStrtR.All_Resp(:,1),1),1) = vpmMockEVStrtR.All_Resp(:,1)
    FRcontEVMock(1:size(vpmMockEVStrtR.rPot(:,1),1),2) = vpmMockEVStrtR.rPot(:,1)
    FRcontEVMock(1:size(vpmMockEVStrtR.rDep(:,1),1),3) = vpmMockEVStrtR.rDep(:,1)
   %Afterinduction
   %l6
   FRAIEVL6=nan(318,3)
    FRAIEVL6(1:size(vpmL6Evok.evokAll(:,1),1),1) = vpmL6Evok.evokAll(:,2)
    FRAIEVL6(1:size(vpmL6Evok.rPot(:,1),1),2) = vpmL6Evok.rPot(:,2)
    FRAIEVL6(1:size(vpmL6Evok.rDep(:,1),1),3) = vpmL6Evok.rDep(:,2)
    %mock
    FRAIEVMock=nan(318,3)
    FRAIEVMock(1:size(vpmMockEVStrtR.All_Resp(:,1),1),1) = vpmMockEVStrtR.All_Resp(:,2)
    FRAIEVMock(1:size(vpmMockEVStrtR.rPot(:,1),1),2) = vpmMockEVStrtR.rPot(:,2)
    FRAIEVMock(1:size(vpmMockEVStrtR.rDep(:,1),1),3) = vpmMockEVStrtR.rDep(:,2)
    %% spontaneous 
    %control l6
    FRcontSPL6=nan(318,3)
    FRcontSPL6(1:size(vpmL6spont.evokAll(:,1),1),1) = vpmL6spont.evokAll(:,1)
    FRcontSPL6(1:size(vpmL6spont.rPot(:,1),1),2) = vpmL6spont.rPot(:,1)
    FRcontSPL6(1:size(vpmL6spont.rDep(:,1),1),3) = vpmL6spont.rDep(:,1)
    %control Mock
    FRcontSPMock=nan(318,3)
    FRcontSPMock(1:size(vpmMockSPstrR.All_Resp(:,1),1),1) = vpmMockSPstrR.All_Resp(:,1)
    FRcontSPMock(1:size(vpmMockSPstrR.rPot(:,1),1),2) = vpmMockSPstrR.rPot(:,1)
    FRcontSPMock(1:size(vpmMockSPstrR.rDep(:,1),1),3) = vpmMockSPstrR.rDep(:,1)
        %afterinduction L6
    FRAISPL6=nan(318,3)
    FRAISPL6(1:size(vpmL6spont.evokAll(:,1),1),1) = vpmL6spont.evokAll(:,1)
    FRAISPL6(1:size(vpmL6spont.rPot(:,1),1),2) = vpmL6spont.rPot(:,1)
    FRAISPL6(1:size(vpmL6spont.rDep(:,1),1),3) = vpmL6spont.rDep(:,1)
    % afterInduction Mock
    FRAISPMock=nan(318,3)
    FRAISPMock(1:size(vpmMockSPstrR.All_Resp(:,1),1),1) = vpmMockSPstrR.All_Resp(:,2)
    FRAISPMock(1:size(vpmMockSPstrR.rPot(:,1),1),2) = vpmMockSPstrR.rPot(:,2)
    FRAISPMock(1:size(vpmMockSPstrR.rDep(:,1),1),3) = vpmMockSPstrR.rDep(:,2)
    
    
    %ALL EVoked and spontaneous Matrix mock&L6
    FREVCAIMOCK_L6=cat(2,FRcontEVL6,FRcontEVMock,FRAIEVL6,FRAIEVMock);
    FRSPCAIMOCK_L6=cat(2,FRcontSPL6,FRcontSPMock,FRAISPL6,FRAIEVMock);
    %names for Evoked groups
    nameEVCL6=cat(2,repmat({'L6contrlEvokALL'},318,1),repmat({'L6contrlEvokPOT'},318,1)...
        ,repmat({'L6contrlEvokDEP'},318,1))
    nameEVAIL6=cat(2,repmat({'L6AfterIEvokALL'},318,1),repmat({'L6AfterIEvokPOT'},318,1)...
        ,repmat({'L6AfterIEvokDEP'},318,1))
    nameEVCMock=cat(2,repmat({'MockContolEvokALL'},318,1),repmat({'MockcontrlEvokPOT'},318,1)...
        ,repmat({'MockcontrlEvokDEP'},318,1))
    nameEVAIMock=cat(2,repmat({'MockAfterIEvokALL'},318,1),repmat({'MockAfterIEvokPOT'},318,1)...
        ,repmat({'MockAfterIEvokDEP'},318,1))
    %concat names evoked
    nameall=cat(2,nameEVCL6,nameEVCMock,nameEVAIL6,nameEVAIMock)
    %names for spontaneous groups
      nameSPCL6=cat(2,repmat({'L6contrlSPontALL'},318,1),repmat({'L6contrlSPpontPOT'},318,1)...
        ,repmat({'L6contrlSPontDEP'},318,1))
      nameSPCMock=cat(2,repmat({'MockContrlSPontALL'},318,1),repmat({'MockcontrlSPontPOT'},318,1)...
        ,repmat({'MockcontrolSPontDEP'},318,1))
     nameSPAIL6=cat(2,repmat({'L6AfterISPontALL'},318,1),repmat({'L6AfterISPontPOT'},318,1)...
        ,repmat({'L6AfterISPontDEP'},318,1))
     nameSPAIMock=cat(2,repmat({'MockAfterISPontALL'},318,1),repmat({'MockAfterISPontPOT'},318,1)...
        ,repmat({'MockAfterISPontDEP'},318,1))
    nameSPall=cat(2,nameSPCL6,nameSPCMock,nameSPAIL6,nameSPAIMock)
    
    FREV_SPCAIMOCK_L6=cat(2,FREVCAIMOCK_L6,FRSPCAIMOCK_L6)
    nameSP_EVall=cat(2,nameall,nameSPall);
    
     [p h stats]=kruskalwallis(FREV_SPCAIMOCK_L6(:),nameSP_EVall(:),'on') %kruskal wallis
     [c m h]=multcompare(stats,'CType','dunn-sidak')
     
    MatSign=nan(6,6)
    for i=1:size(c,1)
        a=c(i,1);b=c(i,2);
        MatSign(a,b)=c(i,6)
    end
     for i=1:size(c,1)
        a=c(i,1);b=c(i,2);
        MatSign(a,b)=c(i,6)
    end
    %%
    %MatSign=nan(24,24);
    for i=1:size(c,1)
        a=c(i,1);b=c(i,2);
        MatSign(b,a)=c(i,6)
    end
    MatSignStep = (MatSign < 0.05) + (MatSign < 0.01) + (MatSign < 0.001);
    range(MatSignStep)
    %% Diff change statistics
    %L6 evoked
    DiffEVAllL6=vpmL6Evok.evokAll(:,2)-vpmL6Evok.evokAll(:,1);
    DiffEVPotL6=vpmL6Evok.rPot(:,2)-vpmL6Evok.rPot(:,1);
    DIffEVDepL6=vpmL6Evok.rDep(:,2)-vpmL6Evok.rDep(:,1);
    %Mock evoked
    DiffEVAllMock=vpmMockEVStrtR.All_Resp(:,2)-vpmMockEVStrtR.All_Resp(:,1);
    DiffEVPotMock=vpmMockEVStrtR.rPot(:,2)-vpmMockEVStrtR.rPot(:,1);
    DIffEVDepMock=vpmMockEVStrtR.rDep(:,2)-vpmMockEVStrtR.rDep(:,1);
    
    %Spontaneous L6
    DiffSPallL6=vpmL6spont.evokAll(:,2)-vpmL6spont.evokAll(:,1);
    DiffSPPotL6=vpmL6spont.rPot(:,2)-vpmL6spont.rPot(:,1);
    DiffSPDepL6=vpmL6spont.rDep(:,2)-vpmL6spont.rDep(:,1);
    
    %spontaneous mock
    
    DiffSPallMock=vpmMockSPstrR.All_Resp(:,2)-vpmMockSPstrR.All_Resp(:,1);
    DiffSPPotMock=vpmMockSPstrR.rPot(:,2)-vpmMockSPstrR.rPot(:,1);
    DiffSPDepMock=vpmMockSPstrR.rDep(:,2)-vpmMockSPstrR.rDep(:,1);
    
    DiffEvSpL6Mock =nan(318,12)
    DiffEvSpL6Mock(1:size(DiffEVAllL6),1)=DiffEVAllL6
    DiffEvSpL6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvSpL6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvSpL6Mock(1:size(DiffEVAllMock),4)=DiffEVAllMock
    DiffEvSpL6Mock(1:size(DiffEVPotMock),5)=DiffEVPotMock
    DiffEvSpL6Mock(1:size(DIffEVDepMock),6)=DIffEVDepMock
    DiffEvSpL6Mock(1:size(DiffSPallL6),7)=DiffSPallL6
    DiffEvSpL6Mock(1:size(DiffSPPotL6),8)=DiffSPPotL6
    DiffEvSpL6Mock(1:size(DiffSPDepL6),9)=DiffSPDepL6
    DiffEvSpL6Mock(1:size(DiffSPallMock),10)=DiffSPallMock
    DiffEvSpL6Mock(1:size(DiffSPPotMock),11)=DiffSPPotMock
    DiffEvSpL6Mock(1:size(DiffSPDepMock),12)=DiffSPDepMock
    %only evoked
    DiffEvL6Mock=nan(318,6)
    DiffEvL6Mock(1:size(DiffEVAllL6),1)=DiffEVAllL6
    DiffEvL6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvL6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvL6Mock(1:size(DiffSPallL6),4)=DiffSPallL6
    DiffEvL6Mock(1:size(DiffSPPotL6),5)=DiffSPPotL6
    DiffEvL6Mock(1:size(DiffSPDepL6),6)=DiffSPDepL6
    %only Spontaneous
    DiffSpL6Mock=nan(318,6)
    DiffSpL6Mock(1:size(DiffSPallL6),1)=DiffSPallL6
    DiffSpL6Mock(1:size(DiffSPPotL6),2)=DiffSPPotL6
    DiffSpL6Mock(1:size(DiffSPDepL6),3)=DiffSPDepL6
    DiffSpL6Mock(1:size(DiffSPallMock),4)=DiffSPallMock
    DiffSpL6Mock(1:size(DiffSPPotMock),5)=DiffSPPotMock
    DiffSpL6Mock(1:size(DiffSPDepMock),6)=DiffSPDepMock
    
    
    
    
       %name diff evoked l6
        nameEVdiffL6=cat(2,repmat({'L6EVokDiffALL'},318,1),repmat({'L6EvokDiffPOT'},318,1)...
        ,repmat({'L6EVokDiffDEP'},318,1))
    %name diff evoked mock
    nameEVdiffMock=cat(2,repmat({'MockEvokDiffALL'},318,1),repmat({'MockEvokDiffPOT'},318,1)...
        ,repmat({'MockEvokDiffDEP'},318,1))
    %name diff sp l6
        nameSPdiffL6=cat(2,repmat({'L6SpontDiffALL'},318,1),repmat({'L6SpontDiffPOT'},318,1)...
        ,repmat({'L6SpontDiffDEP'},318,1))
    %name diff sp mock
       nameSPdiffMock=cat(2,repmat({'MockSPontDiffALL'},318,1),repmat({'MockSpontDiffPOT'},318,1)...
        ,repmat({'MockSpontDiffDEP'},318,1))
    
    nameDiffAll=cat(2,nameEVdiffL6,nameEVdiffMock,nameSPdiffL6,nameSPdiffMock)
    %% DiffEvSPL6Mock all spont and evok comparisons
        MatSignwilkEVSP=nan(12,12);
    for x = 1:size(DiffEvSpL6Mock,2)
        a=DiffEvSpL6Mock(:,x);
        for y = x+1:size(DiffEvSpL6Mock,2)
            b=DiffEvSpL6Mock(:,y);
            fprintf(1, 'x:%d, y:%d\n', x, y)
            MatSignwilkEVSP(y,x)=ranksum(a,b,'tail','right')
        end
    end
    MatSignStepEVSP= (MatSignwilkEVSP<=0.05) + (MatSignwilkEVSP<=0.01) + (MatSignwilkEVSP<=0.001)

NodatEVSP=isnan(MatSignwilkEVSP)| MatSignwilkEVSP>=0.05001
MatSignStepEVSP(NodatEVSP)=NaN

h1=heatmap(MatSignStepEVSP,'YDisplayLabels',nameDiffAll(1,:),'XDisplayLabels',nameDiffAll(1,:),'MissingDataColor','w','GridVisible','off')



testCOlor=[0 0 0;0.2 0.2 0.2; 0.4 0.4 0.4]
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Paired ranksum comparisons L6-VPM vs Mock-VPM')
options=struct(gca)
cb1=options.Colorbar; 
cb1.Ticks=[1.3 2 2.7]
cb1.TickLabels={'*','**','***'}
ax3=options.Axes
ax3.YDir='normal'
       
    
    %%
    MatSignwilkEV=nan(6,6);
    for x = 1:size(DiffEvL6Mock,2)
        a=DiffEvL6Mock(:,x);
        for y = x+1:size(DiffEvL6Mock,2)
            b=DiffEvL6Mock(:,y);
            fprintf(1, 'x:%d, y:%d\n', x, y)
            MatSignwilkEV(y,x)=ranksum(a,b,'tail','right')
        end
    end
 %%   
    MatSignwilkSP=nan(6,6);
    for r = 1:size(DiffSpL6Mock,2)
        c=DiffSpL6Mock(:,r);
        for v = r+1:size(DiffSpL6Mock,2)
            d=DiffSpL6Mock(:,v);
            fprintf(1, 'r:%d, v:%d\n', r, v)
            MatSignwilkSP(v,r)=ranksum(c,d,'tail','right')
        end
    end
        
 

    
MatSignStepEV= (MatSignwilkEV<=0.05) + (MatSignwilkEV<=0.01) + (MatSignwilkEV<=0.001)
Nodat=isnan(MatSignwilkEV)| MatSignwilkEV>=0.05001
MatSignStepEV(Nodat)=NaN

MatSignStepSP= (MatSignwilkSP<=0.05) + (MatSignwilkSP<=0.01) + (MatSignwilkSP<=0.001)
NodatSp=isnan(MatSignwilkSP)| MatSignwilkSP>=0.05001
MatSignStepSP(NodatSp)=NaN

nameDiffAllEv=nameDiffAll(1,[1:6]);
nameDiffAllSP=nameDiffAll(1,[7:12]);

figure
subplot(1,2,1)
h1=heatmap(MatSignStepEV,'YDisplayLabels',nameDiffAllEv,'XDisplayLabels',nameDiffAllEv,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 1 0;0 0.2 0.2; 1 0 0] 
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Evoked Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 

colorbar off
ax3=options.Axes
ax3.YDir='normal'
subplot(1,2,2)
h2=heatmap(MatSignStepSP,'YDisplayLabels',nameDiffAllSP,'XDisplayLabels',nameDiffAllSP,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 0 0;0.2 0.2 0.2; 0.4 0.4 0.4]
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Spontaneous Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 
cb1.Ticks=[0 1 4]
cb1.TickLabels={'*','**','***'}
colorbar off
ax3=options.Axes
ax3.YDir='normal'
%% now with non modulated clusters
%% Diff change statistics
    %L6 evoked
    DiffEVAllL6=vpmL6Evok.evokAll(:,2)-vpmL6Evok.evokAll(:,1);
    DiffEVPotL6=vpmL6Evok.rPot(:,2)-vpmL6Evok.rPot(:,1);
    DIffEVDepL6=vpmL6Evok.rDep(:,2)-vpmL6Evok.rDep(:,1);
    %Mock evoked
    DiffEVAllMock=vpmMockEVStrtR.All_Resp(:,2)-vpmMockEVStrtR.All_Resp(:,1);
    DiffEVPotMock=vpmMockEVStrtR.rPot(:,2)-vpmMockEVStrtR.rPot(:,1);
    DIffEVDepMock=vpmMockEVStrtR.rDep(:,2)-vpmMockEVStrtR.rDep(:,1);
    
    %Spontaneous L6
    DiffSPallL6=vpmL6spont.evokAll(:,2)-vpmL6spont.evokAll(:,1);
    DiffSPPotL6=vpmL6spont.rPot(:,2)-vpmL6spont.rPot(:,1);
    DiffSPDepL6=vpmL6spont.rDep(:,2)-vpmL6spont.rDep(:,1);
    
    %spontaneous mock
    
    DiffSPallMock=vpmMockSPstrR.All_Resp(:,2)-vpmMockSPstrR.All_Resp(:,1);
    DiffSPPotMock=vpmMockSPstrR.rPot(:,2)-vpmMockSPstrR.rPot(:,1);
    DiffSPDepMock=vpmMockSPstrR.rDep(:,2)-vpmMockSPstrR.rDep(:,1);
    
    DiffEvSpL6Mock =nan(318,12)
    DiffEvSpL6Mock(1:size(DiffEVAllL6),1)=DiffEVAllL6
    DiffEvSpL6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvSpL6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvSpL6Mock(1:size(DiffEVAllMock),4)=DiffEVAllMock
    DiffEvSpL6Mock(1:size(DiffEVPotMock),5)=DiffEVPotMock
    DiffEvSpL6Mock(1:size(DIffEVDepMock),6)=DIffEVDepMock
    DiffEvSpL6Mock(1:size(DiffSPallL6),7)=DiffSPallL6
    DiffEvSpL6Mock(1:size(DiffSPPotL6),8)=DiffSPPotL6
    DiffEvSpL6Mock(1:size(DiffSPDepL6),9)=DiffSPDepL6
    DiffEvSpL6Mock(1:size(DiffSPallMock),10)=DiffSPallMock
    DiffEvSpL6Mock(1:size(DiffSPPotMock),11)=DiffSPPotMock
    DiffEvSpL6Mock(1:size(DiffSPDepMock),12)=DiffSPDepMock
    %only evoked
    DiffEvL6Mock=nan(318,6)
    DiffEvL6Mock(1:size(DiffEVAllL6),1)=DiffEVAllL6
    DiffEvL6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvL6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvL6Mock(1:size(DiffEVAllMock),4)=DiffEVAllMock
    DiffEvL6Mock(1:size(DiffEVPotMock),5)=DiffEVPotMock
    DiffEvL6Mock(1:size(DIffEVDepMock),6)=DIffEVDepMock
    %only Spontaneous
    DiffSpL6Mock=nan(318,6)
    DiffSpL6Mock(1:size(DiffSPallL6),1)=DiffSPallL6
    DiffSpL6Mock(1:size(DiffSPPotL6),2)=DiffSPPotL6
    DiffSpL6Mock(1:size(DiffSPDepL6),3)=DiffSPDepL6
    DiffSpL6Mock(1:size(DiffSPallMock),4)=DiffSPallMock
    DiffSpL6Mock(1:size(DiffSPPotMock),5)=DiffSPPotMock
    DiffSpL6Mock(1:size(DiffSPDepMock),6)=DiffSPDepMock
    
    
    
    
       %name diff evoked l6
        nameEVdiffL6=cat(2,repmat({'L6EVokDiffALL'},318,1),repmat({'L6EvokDiffPOT'},318,1)...
        ,repmat({'L6EVokDiffDEP'},318,1))
    %name diff evoked mock
    nameEVdiffMock=cat(2,repmat({'MockEvokDiffALL'},318,1),repmat({'MockEvokDiffPOT'},318,1)...
        ,repmat({'MockEvokDiffDEP'},318,1))
    %name diff sp l6
        nameSPdiffL6=cat(2,repmat({'L6SpontDiffALL'},318,1),repmat({'L6SpontDiffPOT'},318,1)...
        ,repmat({'L6SpontDiffDEP'},318,1))
    %name diff sp mock
       nameSPdiffMock=cat(2,repmat({'MockSPontDiffALL'},318,1),repmat({'MockSpontDiffPOT'},318,1)...
        ,repmat({'MockSpontDiffDEP'},318,1))
    
    nameDiffAll=cat(2,nameEVdiffL6,nameEVdiffMock,nameSPdiffL6,nameSPdiffMock)
    %% DiffEvSPL6Mock all spont and evok comparisons
        MatSignwilkEVSP=nan(12,12);
    for x = 1:size(DiffEvSpL6Mock,2)
        a=DiffEvSpL6Mock(:,x);
        for y = x+1:size(DiffEvSpL6Mock,2)
            b=DiffEvSpL6Mock(:,y);
            fprintf(1, 'x:%d, y:%d\n', x, y)
            MatSignwilkEVSP(y,x)=ranksum(a,b,'tail','right')
        end
    end
    MatSignStepEVSP= (MatSignwilkEVSP<=0.05) + (MatSignwilkEVSP<=0.01) + (MatSignwilkEVSP<=0.001)

NodatEVSP=isnan(MatSignwilkEVSP)| MatSignwilkEVSP>=0.05001
MatSignStepEVSP(NodatEVSP)=NaN

h1=heatmap(MatSignStepEVSP,'YDisplayLabels',nameDiffAll(1,:),'XDisplayLabels',nameDiffAll(1,:),'MissingDataColor','w','GridVisible','off')



testCOlor=[0 0 0;0.2 0.2 0.2; 0.4 0.4 0.4]
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Paired ranksum comparisons L6-VPM vs Mock-VPM')
options=struct(gca)
cb1=options.Colorbar; 
cb1.Ticks=[1.3 2 2.7]
cb1.TickLabels={'*','**','***'}
ax3=options.Axes
ax3.YDir='normal'
       
    
    %%
    MatSignwilkEV=nan(6,6);
    for x = 1:size(DiffEvL6Mock,2)
        a=DiffEvL6Mock(:,x);
        for y = x+1:size(DiffEvL6Mock,2)
            b=DiffEvL6Mock(:,y);
            fprintf(1, 'x:%d, y:%d\n', x, y)
            MatSignwilkEV(y,x)=ranksum(a,b,'tail','right')
        end
    end
 %%   
    MatSignwilkSP=nan(6,6);
    for r = 1:size(DiffSpL6Mock,2)
        c=DiffSpL6Mock(:,r);
        for v = r+1:size(DiffSpL6Mock,2)
            d=DiffSpL6Mock(:,v);
            fprintf(1, 'r:%d, v:%d\n', r, v)
            MatSignwilkSP(v,r)=ranksum(c,d,'tail','right')
        end
    end
        
 

    
MatSignStepEV= (MatSignwilkEV<=0.05) + (MatSignwilkEV<=0.01) + (MatSignwilkEV<=0.001)
Nodat=isnan(MatSignwilkEV)| MatSignwilkEV>=0.05001
MatSignStepEV(Nodat)=NaN

MatSignStepSP= (MatSignwilkSP<=0.05) + (MatSignwilkSP<=0.01) + (MatSignwilkSP<=0.001)
NodatSp=isnan(MatSignwilkSP)| MatSignwilkSP>=0.05001
MatSignStepSP(NodatSp)=NaN

nameDiffAllEv=nameDiffAll(1,[1:6]);
nameDiffAllSP=nameDiffAll(1,[7:12]);

figure
subplot(1,2,1)
h1=heatmap(MatSignStepEV,'YDisplayLabels',nameDiffAllEv,'XDisplayLabels',nameDiffAllEv,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 1 0;0 0.2 0.2; 1 0 0] 
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Evoked Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 

colorbar off
ax3=options.Axes
ax3.YDir='normal'
subplot(1,2,2)
h2=heatmap(MatSignStepSP,'YDisplayLabels',nameDiffAllSP,'XDisplayLabels',nameDiffAllSP,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 0 0;0.2 0.2 0.2; 0.4 0.4 0.4]
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Spontaneous Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 
cb1.Ticks=[0 1 4]
cb1.TickLabels={'*','**','***'}
colorbar off
ax3=options.Axes
ax3.YDir='normal'
%% now with non modulated clusters
%Diff change statistics
    %L6 evoked
    DiffEVNonML6=vpmL6EVStrtR2.NonMOd(:,2)-vpmL6EVStrtR2.NonMOd(:,1);
    DiffEVPotL6=vpmL6Evok.rPot(:,2)-vpmL6Evok.rPot(:,1);
    DIffEVDepL6=vpmL6Evok.rDep(:,2)-vpmL6Evok.rDep(:,1);
    %Mock evoked
    DiffEVNonMMock=vpmMockEVStrtR2.NonMod(:,2)-vpmMockEVStrtR2.NonMod(:,1);
    DiffEVPotMock=vpmMockEVStrtR.rPot(:,2)-vpmMockEVStrtR.rPot(:,1);
    DIffEVDepMock=vpmMockEVStrtR.rDep(:,2)-vpmMockEVStrtR.rDep(:,1);
    
    %Spontaneous L6
    DiffSPNonML6=vpmL6SPstrR2.NonMod(:,2)-vpmL6SPstrR2.NonMod(:,1);
    DiffSPPotL6=vpmL6spont.rPot(:,2)-vpmL6spont.rPot(:,1);
    DiffSPDepL6=vpmL6spont.rDep(:,2)-vpmL6spont.rDep(:,1);
    
    %spontaneous mock
    
    DiffSPNonMMock=vpmMockSPstrR2.NonMod(:,2)-vpmMockSPstrR2.NonMod(:,1);
    DiffSPPotMock=vpmMockSPstrR.rPot(:,2)-vpmMockSPstrR.rPot(:,1);
    DiffSPDepMock=vpmMockSPstrR.rDep(:,2)-vpmMockSPstrR.rDep(:,1);
    
    DiffEvSpWnML6Mock =nan(318,12)
    DiffEvSpWnML6Mock(1:size(DiffEVNonML6),1)=DiffEVNonML6
    DiffEvSpWnML6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvSpWnML6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvSpWnML6Mock(1:size(DiffEVNonMMock),4)=DiffEVNonMMock
    DiffEvSpWnML6Mock(1:size(DiffEVPotMock),5)=DiffEVPotMock
    DiffEvSpWnML6Mock(1:size(DIffEVDepMock),6)=DIffEVDepMock
    DiffEvSpWnML6Mock(1:size(DiffSPNonML6),7)=DiffSPNonML6
    DiffEvSpWnML6Mock(1:size(DiffSPPotL6),8)=DiffSPPotL6
    DiffEvSpWnML6Mock(1:size(DiffSPDepL6),9)=DiffSPDepL6
    DiffEvSpWnML6Mock(1:size(DiffSPNonMMock),10)=DiffSPNonMMock
    DiffEvSpWnML6Mock(1:size(DiffSPPotMock),11)=DiffSPPotMock
    DiffEvSpWnML6Mock(1:size(DiffSPDepMock),12)=DiffSPDepMock
    %only evoked
    DiffEvWnML6Mock=nan(318,6)
    DiffEvWnML6Mock(1:size(DiffEVNonML6),1)=DiffEVNonML6
    DiffEvWnML6Mock(1:size(DiffEVPotL6),2)=DiffEVPotL6
    DiffEvWnML6Mock(1:size(DIffEVDepL6),3)=DIffEVDepL6
    DiffEvWnML6Mock(1:size(DiffEVNonMMock),4)=DiffEVNonMMock
    DiffEvWnML6Mock(1:size(DiffEVPotMock),5)=DiffEVPotMock
    DiffEvWnML6Mock(1:size(DIffEVDepMock),6)=DIffEVDepMock
    %only Spontaneous
    DiffSpL6Mock=nan(318,6)
    DiffSpL6Mock(1:size(DiffSPallL6),1)=DiffSPallL6
    DiffSpL6Mock(1:size(DiffSPPotL6),2)=DiffSPPotL6
    DiffSpL6Mock(1:size(DiffSPDepL6),3)=DiffSPDepL6
    DiffSpL6Mock(1:size(DiffSPNonMMock),4)=DiffSPNonMMock
    DiffSpL6Mock(1:size(DiffSPPotMock),5)=DiffSPPotMock
    DiffSpL6Mock(1:size(DiffSPDepMock),6)=DiffSPDepMock
    
    
    
    
       %name diff evoked l6
        nameEVdiffL6=cat(2,repmat({'L6EVokDiffNonM'},318,1),repmat({'L6EvokDiffPOT'},318,1)...
        ,repmat({'L6EVokDiffDEP'},318,1))
    %name diff evoked mock
    nameEVdiffMock=cat(2,repmat({'MockEvokDiffNonM'},318,1),repmat({'MockEvokDiffPOT'},318,1)...
        ,repmat({'MockEvokDiffDEP'},318,1))
    %name diff sp l6
        nameSPdiffL6=cat(2,repmat({'L6SpontDiffNonM'},318,1),repmat({'L6SpontDiffPOT'},318,1)...
        ,repmat({'L6SpontDiffDEP'},318,1))
    %name diff sp mock
       nameSPdiffMock=cat(2,repmat({'MockSPontDiffNonM'},318,1),repmat({'MockSpontDiffPOT'},318,1)...
        ,repmat({'MockSpontDiffDEP'},318,1))
    
    nameDiffAll=cat(2,nameEVdiffL6,nameEVdiffMock,nameSPdiffL6,nameSPdiffMock)


  MatSignwilkNonMEV=nan(6,6);
    for x = 1:size(DiffEvL6Mock,2)
        a=DiffEvL6Mock(:,x);
        for y = x+1:size(DiffEvL6Mock,2)
            b=DiffEvL6Mock(:,y);
            fprintf(1, 'x:%d, y:%d\n', x, y)
            MatSignwilkNonMEV(y,x)=ranksum(a,b)
        end
    end

   
    MatSignwilkNonMSP=nan(6,6);
    for r = 1:size(DiffSpL6Mock,2)
        c=DiffSpL6Mock(:,r);
        for v = r+1:size(DiffSpL6Mock,2)
            d=DiffSpL6Mock(:,v);
            fprintf(1, 'r:%d, v:%d\n', r, v)
            MatSignwilkNonMSP(v,r)=ranksum(c,d)
        end
    end
        
     
MatSignStepNonMEV= (MatSignwilkNonMEV<=0.05) + (MatSignwilkNonMEV<=0.01) + (MatSignwilkNonMEV<=0.001)
NodatNonM=isnan(MatSignwilkNonMEV)%| MatSignwilkNonMEV>=0.05001

MatSignStepNonMEV(NodatNonM)=NaN

MatSignStepNonMSP= (MatSignwilkNonMSP<=0.05) + (MatSignwilkNonMSP<=0.01) + (MatSignwilkNonMSP<=0.001)
NodatNonMSp=isnan(MatSignwilkNonMSP)%| MatSignwilkNonMSP>=0.05001


MatSignStepNonMSP(NodatNonMSp)=NaN

nameDiffAllEv=nameDiffAll(1,[1:6]);
nameDiffAllSP=nameDiffAll(1,[7:12]);

figure
subplot(1,2,1)
h1=heatmap(MatSignStepNonMEV,'YDisplayLabels',nameDiffAllEv,'XDisplayLabels',nameDiffAllEv,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 1 0;0 0.2 0.2; 1 0 0] 
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Evoked Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 

colorbar off
ax3=options.Axes
ax3.YDir='normal'
subplot(1,2,2)
h2=heatmap(MatSignStepNonMSP,'YDisplayLabels',nameDiffAllSP,'XDisplayLabels',nameDiffAllSP,...
    'MissingDataColor','w','GridVisible','off')
testCOlor=[0 0 0;0.2 0.2 0.2; 0.4 0.4 0.4]
colormap(testCOlor)
%ax=axes;
%ax.Color='none'
%plot(
title('Spontaneous Paired ranksum comparisons L6 vs Mock')
options=struct(gca)
cb1=options.Colorbar; 
cb1.Ticks=[0 1 4]
cb1.TickLabels={'*','**','***'}
colorbar off
ax3=options.Axes
ax3.YDir='normal'
    
    
    
    