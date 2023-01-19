%% Violing plot for structure Mock vs real L6 induction difference AfterInduction - control
DiffMockEVPot=vpmMockEVStrtR.rPot(:,2)-vpmMockEVStrtR.rPot(:,1)
DiffL6EVPot=vpmL6Evok.rPot(:,2)-vpmL6Evok.rPot(:,1)

h2=subplot(2,1,2)
[yCtrE xCtrE]=ksdensity(DiffMockEVPot,'Bandwidth',7)
patch(yCtrE*-1,xCtrE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIE xAIE]=ksdensity(DiffL6EVPot,'Bandwidth',7)
patch(yAIE,xAIE,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)-0.0025,DiffMockEVPot,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(68,1)+0.0025,DiffL6EVPot,'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.0025,median(DiffMockEVPot),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0025,median(DiffMockEVPot),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.0025,median(DiffL6EVPot),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0025,median(DiffL6EVPot),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
%deprressed
DiffMockEVDep=vpmMockEVStrtR.rDep(:,2)-vpmMockEVStrtR.rDep(:,1)
DiffL6EVDep=vpmL6Evok.rDep(:,2)-vpmL6Evok.rDep(:,1)

[yCtrDepE xCtrDepE]=ksdensity(DiffMockEVDep,'Bandwidth',7)
patch((yCtrDepE-0.075)*-1,xCtrDepE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDepE xAIDepE]=ksdensity(DiffL6EVDep,'Bandwidth',7)
patch((yAIDepE+0.075),xAIDepE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(38,1)+0.0725,DiffMockEVDep,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.0775,DiffL6EVDep,'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.0725,median(DiffMockEVDep),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0725,median(DiffMockEVDep),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.0775,median(DiffL6EVDep),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0775,median(DiffL6EVDep),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)

%all responsive
DiffMockEVAll=vpmMockEVStrtR.All_Resp(:,2)-vpmMockEVStrtR.All_Resp(:,1)
DiffL6EVAll=vpmL6Evok.evokAll(:,2)-vpmL6Evok.evokAll(:,1)


[yCtrAllE xCtrAllE]=ksdensity(DiffMockEVAll,'Bandwidth',7)
patch((yCtrAllE+0.075)*-1,xCtrAllE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAllE xAIAllE]=ksdensity(DiffL6EVAll,'Bandwidth',7)

patch((yAIAllE-0.075),xAIAllE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7) % evoked
hold on
plot(zeros(318,1)-0.0725,DiffL6EVAll,'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(251,1)-0.0775,DiffMockEVAll,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.0725,median(DiffL6EVAll),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0725,median(DiffL6EVAll),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(-0.0775,median(DiffMockEVAll),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0775,median(DiffMockEVAll),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
ylim([-160 160])
xlim([-0.125 0.125])

set(gca,'XTick',[-0.07 0.0 0.07])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[-120 -80 -40  0 40 80 120])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
%titl=title ('      Evoked VPM firing rate  After Induction: \newlineDifference (AI-Control) L6 vs Sham stimulation')
%set(titl,'horizontalAlignment','Center')
xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
%annotation('line',[0.51 0.76],[0.85 0.85],'LineWidth',1);
%text([0.025],[255],{'Significantly Modulated'},'HorizontalAlignment','center')
%legend('Control','After-Induction','Location','northwest'); box off; hold off
legend('Sham Induction','L6 Induction'); box off; 
set(legend,'position',[0.735 0.46 0.21 0.075]); legend boxon

%% Violing plot for structure Mock vs real L6 induction difference AfterInduction - control Spontaneous
DiffMockSPPot=vpmMockSPstrR.rPot(:,2)-vpmMockSPstrR.rPot(:,1)
DiffL6SPPot=vpmL6spont.rPot(:,2)-vpmL6spont.rPot(:,1)

h1=subplot(2,1,1)
h1.XColor=[1 1 1]
[yCtrE xCtrE]=ksdensity(DiffMockSPPot,'Bandwidth',1.5)
patch(yCtrE*-1,xCtrE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIE xAIE]=ksdensity(DiffL6SPPot,'Bandwidth',1.5)
patch(yAIE,xAIE,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)-0.01,DiffMockSPPot,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(68,1)+0.01,DiffL6SPPot,'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.0125,median(DiffMockSPPot),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0125,median(DiffMockSPPot),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.0125,median(DiffL6SPPot),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0125,median(DiffL6SPPot),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
%deprressed
DiffMockSPDep=vpmMockSPstrR.rDep(:,2)-vpmMockSPstrR.rDep(:,1)
DiffL6SPDep=vpmL6spont.rDep(:,2)-vpmL6spont.rDep(:,1)

[yCtrDepE xCtrDepE]=ksdensity(DiffMockSPDep,'Bandwidth',1.5)
patch((yCtrDepE-0.375)*-1,xCtrDepE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDepE xAIDepE]=ksdensity(DiffL6SPDep,'Bandwidth',1.5)
patch((yAIDepE+0.375),xAIDepE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(38,1)+0.365,DiffMockSPDep,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.385,DiffL6SPDep,'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.365,median(DiffMockSPDep),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.365,median(DiffMockSPDep),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.385,median(DiffL6SPDep),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.385,median(DiffL6SPDep),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)

%all responsive
DiffMockSPAll=vpmMockSPstrR.All_Resp(:,2)-vpmMockSPstrR.All_Resp(:,1)
DiffL6SPAll=vpmL6spont.evokAll(:,2)-vpmL6spont.evokAll(:,1)


[yCtrAllE xCtrAllE]=ksdensity(DiffMockSPAll,'Bandwidth',7)
patch((yCtrAllE+0.375)*-1,xCtrAllE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAllE xAIAllE]=ksdensity(DiffL6SPAll,'Bandwidth',7)

patch((yAIAllE-0.375),xAIAllE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7) % evoked
hold on
plot(zeros(318,1)-0.365,DiffL6SPAll,'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(251,1)-0.385,DiffMockEVAll,'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.365,median(DiffL6EVAll),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.365,median(DiffL6EVAll),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(-0.385,median(DiffMockEVAll),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.385,median(DiffMockEVAll),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
ylim([-50 50])
xlim([-0.625 0.625])
set(gca,'XTick',[-0.375 0.0 0.375])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[-50 -25 0 25 50])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
%titl=title ('      Evoked VPM firing rate  After Induction: \newlineDifference (AI-Control) L6 vs Sham stimulation')
%set(titl,'horizontalAlignment','Center')
%xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
%annotation('line',[0.51 0.76],[0.85 0.85],'LineWidth',1);
%text([0.025],[255],{'Significantly Modulated'},'HorizontalAlignment','center')
%legend('Control','After-Induction','Location','northwest'); box off; hold off
legend('Sham Induction','L6 Induction'); box off; 
set(legend,'position',[0.735 0.46 0.21 0.075]); legend boxon
set(gcf,'color','w');
hold off
%%
ylim([-100 70])
annotation('line',[0.505 0.53],[0.81 0.81],'LineWidth',1);
annotation('line',[0.735 0.76],[0.81 0.81],'LineWidth',1);
annotation('line',[0.2725 0.2975],[0.81 0.81],'LineWidth',1);
text([0],[50],{'ns'},'HorizontalAlignment','center','fontsize',14)
text([-0.375],[50],{'***'},'HorizontalAlignment','center','fontsize',14)
text([0.375],[50],{'ns'},'HorizontalAlignment','center','fontsize',14)