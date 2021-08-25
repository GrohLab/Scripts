
[yCtrE xCtrE]=ksdensity(evFr(wruIdx&posMod&signMod,1),'Bandwidth',7)
patch(yCtrE*-1,xCtrE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIE xAIE]=ksdensity(evFr(wruIdx&posMod&signMod,2),'Bandwidth',7)
patch(yAIE,xAIE,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(68,1)-0.0025,evFr(wruIdx&posMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(68,1)+0.0025,evFr(wruIdx&posMod&signMod,2),'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.0025,median(evFr(wruIdx&posMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.0025,median(evFr(wruIdx&posMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(+0.0025,median(evFr(wruIdx&posMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.0025,median(evFr(wruIdx&posMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)

%deprressed

[yCtrDepE xCtrDepE]=ksdensity(evFr(wruIdx&negMod&signMod,1),'Bandwidth',7)
patch((yCtrDepE-0.05)*-1,xCtrDepE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDepE xAIDepE]=ksdensity(evFr(wruIdx&negMod&signMod,2),'Bandwidth',7)
patch((yAIDepE+0.05),xAIDepE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)+0.0475,evFr(wruIdx&negMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.0525,evFr(wruIdx&negMod&signMod,2),'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.0475,median(evFr(wruIdx&negMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.0475,median(evFr(wruIdx&negMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(+0.0525,median(evFr(wruIdx&negMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.0525,median(evFr(wruIdx&negMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)

%all responsive

[yCtrAllE xCtrAllE]=ksdensity(evFr(wruIdx,1),'Bandwidth',7)
patch((yCtrAllE+0.05)*-1,xCtrAllE,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAllE xAIAllE]=ksdensity(evFr(wruIdx,2),'Bandwidth',7)
patch((yAIAllE-0.05),xAIAllE+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(318,1)-0.0475,evFr(wruIdx,2),'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(318,1)-0.0525,evFr(wruIdx,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.0475,median(evFr(wruIdx,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.0475,median(evFr(wruIdx,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(-0.0525,median(evFr(wruIdx,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.0525,median(evFr(wruIdx,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
ylim([-50 300])
xlim([-0.08 0.08])

set(gca,'XTick',[-0.05 0.0 0.05])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[0,50,100,150 200])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
title 'Evoked firing rate for VPM clusters'
xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
annotation('line',[0.51 0.76],[0.8 0.8],'LineWidth',1);
text([0.024],[255],{'Significantly Modulated'},'HorizontalAlignment','center')
legend('Control','After-Induction','Location','northwest'); box off; hold off