plot(zeros(68,1)-0.05,pfr(wruIdx&posMod&signMod,1),'ksquare','MarkerSize',12,'MarkerFaceColor',[0 0.6 1])
hold on
plot(zeros(68,1)+0.05,pfr(wruIdx&posMod&signMod,2),'k^','MarkerSize',10,'MarkerFaceColor',[1 0 0])

plot(zeros(33,1)+0.45,pfr(wruIdx&negMod&signMod,1),'ksquare','MarkerSize',12,'MarkerFaceColor',[0 0.6 1])
xlim([-1 1])
plot(zeros(33,1)+0.55,pfr(wruIdx&negMod&signMod,2),'kv','MarkerSize',10,'MarkerFaceColor',[1 0 0])


plot(zeros(318,1)-0.55,pfr(wruIdx,1),'ko','MarkerSize',10,'MarkerFaceColor',[0 0.6 1])
ylim([-5 40])
plot(zeros(318,1)-0.45,pfr(wruIdx,2),'ko','MarkerSize',10,'MarkerFaceColor',[1 0 0])

%plot(zeros(217,1)+-1.1,pfr(wruIdx&signMod==0,1),'ksquare','MarkerSize',10,'MarkerFaceColor',[0 0.6 1]); hold on
%plot(zeros(217,1)+-1,pfr(wruIdx&signMod==0,2),'ksquare','MarkerSize',10,'MarkerFaceColor',[1 0 0])

%plot(zeros(1928,1)+-1.6,pfr(:,1),'ksquare','MarkerSize',10,'MarkerFaceColor',[0 0.6 1])
%plot(zeros(1928,1)+-1.5,pfr(:,2),'ksquare','MarkerSize',10,'MarkerFaceColor',[1 0 0])

set(gca,'XTick',[-0.5:0.5:0.5])
set(gca,'XTickLabel',["Responsive All", "Potentiated","Depressed"])
title 'spontaneous firing rate for VPM clusters'
xlabel("VPM Clusters")
ylabel("Firing Rate [Hz]")
annotation('line',[0.5 0.75],[0.8 0.8],'LineWidth',1);
text([-0.05],[35],'significantly Modulated','HorizontalAlignment','left')
legend('Control','After-Induction','Location','northwest'); box off; hold off