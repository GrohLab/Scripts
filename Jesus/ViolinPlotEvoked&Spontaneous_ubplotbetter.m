posMod=MIevok>=0.001;negMod=MIevok<=-0.001
h1=subplot(2,1,2)
[yCtrE xCtrE]=ksdensity(evFr(wruIdx&posMod&signMod,1),'Bandwidth',7)
patch(yCtrE*-1,xCtrE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIE xAIE]=ksdensity(evFr(wruIdx&posMod&signMod,2),'Bandwidth',7)
patch(yAIE,xAIE,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)-0.0025,evFr(wruIdx&posMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.0025,evFr(wruIdx&posMod&signMod,2),'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.0025,median(evFr(wruIdx&posMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0025,median(evFr(wruIdx&posMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.0025,median(evFr(wruIdx&posMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0025,median(evFr(wruIdx&posMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)


%deprressed

[yCtrDepE xCtrDepE]=ksdensity(evFr(wruIdx&negMod&signMod,1),'Bandwidth',7)
patch((yCtrDepE-0.05)*-1,xCtrDepE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDepE xAIDepE]=ksdensity(evFr(wruIdx&negMod&signMod,2),'Bandwidth',7)
patch((yAIDepE+0.05),xAIDepE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(38,1)+0.0475,evFr(wruIdx&negMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(38,1)+0.0525,evFr(wruIdx&negMod&signMod,2),'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.0475,median(evFr(wruIdx&negMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0475,median(evFr(wruIdx&negMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.0525,median(evFr(wruIdx&negMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.0525,median(evFr(wruIdx&negMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)

%all responsive

[yCtrAllE xCtrAllE]=ksdensity(evFr(wruIdx,1),'Bandwidth',7)
patch((yCtrAllE+0.05)*-1,xCtrAllE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAllE xAIAllE]=ksdensity(evFr(wruIdx,2),'Bandwidth',7)

patch((yAIAllE-0.05),xAIAllE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7) % evoked
hold on
plot(zeros(251,1)-0.0475,evFr(wruIdx,2),'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(251,1)-0.0525,evFr(wruIdx,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.0475,median(evFr(wruIdx,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0475,median(evFr(wruIdx,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(-0.0525,median(evFr(wruIdx,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.0525,median(evFr(wruIdx,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
ylim([-50 300])
xlim([-0.08 0.08])

set(gca,'XTick',[-0.05 0.0 0.05])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[0 40,80,120,160 200])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
title 'Evoked firing rate for VPM clusters'
xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
%annotation('line',[0.51 0.76],[0.85 0.85],'LineWidth',1);
%text([0.025],[255],{'Significantly Modulated'},'HorizontalAlignment','center')
%legend('Control','After-Induction','Location','northwest'); box off; hold off
legend('Control','After-Induction'); box off; 
set(legend,'position',[0.735 0.46 0.21 0.075]); legend boxon
%% spontaneous 
h2=subplot(2,1,1)
h2.XColor=[1 1 1]

%potentiated

[yCtr xCtr]=ksdensity(pfr(wruIdx&posMod&signMod,1),'Bandwidth',0.7)
patch(yCtr*-1,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAI xAI]=ksdensity(pfr(wruIdx&posMod&signMod,2),'Bandwidth',0.7)
patch(yAI,xAI,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)-0.025,pfr(wruIdx&posMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.025,pfr(wruIdx&posMod&signMod,2),'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.025,median(pfr(wruIdx&posMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.025,median(pfr(wruIdx&posMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.025,median(pfr(wruIdx&posMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.025,median(pfr(wruIdx&posMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)

%deprressed

[yCtrDep xCtrDep]=ksdensity(pfr(wruIdx&negMod&signMod,1),'Bandwidth',1)
patch((yCtrDep-0.5)*-1,xCtrDep,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDep xAIDep]=ksdensity(pfr(wruIdx&negMod&signMod,2),'Bandwidth',1)
patch((yAIDep+0.5),xAIDep+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(38,1)+0.475,pfr(wruIdx&negMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(38,1)+0.525,pfr(wruIdx&negMod&signMod,2),'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.475,median(pfr(wruIdx&negMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.475,median(pfr(wruIdx&negMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(+0.525,median(pfr(wruIdx&negMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(+0.525,median(pfr(wruIdx&negMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)

%all responsive

[yCtrAll xCtrAll]=ksdensity(pfr(wruIdx,1),'Bandwidth',1)
patch((yCtrAll+0.5)*-1,xCtrAll,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAll xAIAll]=ksdensity(pfr(wruIdx,2),'Bandwidth',1)
patch((yAIAll-0.5),xAIAll+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(251,1)-0.475,pfr(wruIdx,2),'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(8,1)-0.525,pfr(wruIdx,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.475,median(pfr(wruIdx,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.475,median(pfr(wruIdx,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
plot(-0.525,median(pfr(wruIdx,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',6)
plot(-0.525,median(pfr(wruIdx,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',6)
ylim([-5 40])
title 'VPM Responsive Clusters Distribution'
set(gca,'XTick',[-0.5:0.5:0.5])
%set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[0,10,20,30])
set(gca,'XTickLabel',[])
title 'Spontaneous firing rate for VPM clusters'
%xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
%annotation('line',[0.51 0.76],[0.385 0.385],'LineWidth',1);
annotation('line',[0.51 0.76],[0.86 0.86],'LineWidth',1);
text([0.24],[35],{'Significantly Modulated'},'HorizontalAlignment','center')
set(gcf,'color','w');
hold off