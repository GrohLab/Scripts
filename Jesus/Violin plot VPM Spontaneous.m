%subplot(2,1,1)
%potentiated

[yCtr xCtr]=ksdensity(pfr(wruIdx&posMod&signMod,1),'Bandwidth',0.7)
patch(yCtr*-1,xCtr,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAI xAI]=ksdensity(pfr(wruIdx&posMod&signMod,2),'Bandwidth',0.7)
patch(yAI,xAI,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(68,1)-0.025,pfr(wruIdx&posMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(68,1)+0.025,pfr(wruIdx&posMod&signMod,2),'k^','MarkerFaceColor',[0.7 0 0])
plot(-0.025,median(pfr(wruIdx&posMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.025,median(pfr(wruIdx&posMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(+0.025,median(pfr(wruIdx&posMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.025,median(pfr(wruIdx&posMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)

%deprressed

[yCtrDep xCtrDep]=ksdensity(pfr(wruIdx&negMod&signMod,1),'Bandwidth',1)
patch((yCtrDep-0.5)*-1,xCtrDep,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIDep xAIDep]=ksdensity(pfr(wruIdx&negMod&signMod,2),'Bandwidth',1)
patch((yAIDep+0.5),xAIDep+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(33,1)+0.475,pfr(wruIdx&negMod&signMod,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(zeros(33,1)+0.525,pfr(wruIdx&negMod&signMod,2),'kv','MarkerFaceColor',[0.7 0 0])
plot(+0.475,median(pfr(wruIdx&negMod&signMod,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.475,median(pfr(wruIdx&negMod&signMod,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(+0.525,median(pfr(wruIdx&negMod&signMod,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(+0.525,median(pfr(wruIdx&negMod&signMod,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)

%all responsive

[yCtrAll xCtrAll]=ksdensity(pfr(wruIdx,1),'Bandwidth',1)
patch((yCtrAll+0.5)*-1,xCtrAll,[0 0 1],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
[yAIAll xAIAll]=ksdensity(pfr(wruIdx,2),'Bandwidth',1)
patch((yAIAll-0.5),xAIAll+1,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
hold on
plot(zeros(318,1)-0.475,pfr(wruIdx,2),'ksquare','MarkerFaceColor',[0.7 0 0])
plot(zeros(318,1)-0.525,pfr(wruIdx,1),'ksquare','MarkerFaceColor',[0 0 0.7])
plot(-0.475,median(pfr(wruIdx,2)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.475,median(pfr(wruIdx,2)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
plot(-0.525,median(pfr(wruIdx,1)),'w*','MarkerFaceColor',[0 0 0],'MarkerSize',8)
plot(-0.525,median(pfr(wruIdx,1)),'ko','MarkerFaceColor',[1 1 1],'MarkerSize',8)
ylim([-5 40])
title 'VPM Responsive Clusters Distribution'
set(gca,'XTick',[-0.5:0.5:0.5])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
set(gca,'yTick',[0,10,20,30])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
title 'Spontaneous firing rate for VPM clusters'
xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
annotation('line',[0.51 0.76],[0.8 0.8],'LineWidth',1);
text([0.24],[35],{'Significantly Modulated'},'HorizontalAlignment','center')
legend('Control','After-Induction','Location','northwest'); box off; hold off
%% pie 

axes('Position',[0.8 0.8 0.15 0.15])
pie([68 33])
colormap([0.25 0.25 0.25; 0.5 0.5 0.5])
box on




%% C=colormap(hot(100))
patch(yCtr*-1,xCtr,[1 0 0],'EdgeAlpha',0.1,'FaceAlpha',0.7)
fill(yAI,xAI,C)
colormap(hot(100))
colorbar