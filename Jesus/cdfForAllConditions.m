x
%CdfOff response

figure(71);

cdfplot(allControlOff); hold on;
cdfplot(allL200Off); hold on;
cdfplot(allL100Off); hold on;
cdfplot(allL50Off); hold on
cdfplot(allL10Off); hold on;
cdfplot(allL1Off); legend("Control","L200","L100","L50","L10","L1")
xlabel('time(ms)');ylabel('probability')
title('Cumulative Distribution all conditions Off')
saveas(figure(71),[pwd '/New figures Jesus/all/cdfOffAll.emf']);
saveas(figure(71),[pwd '/New figures Jesus/all/cdfOffAll.pdf']);

%Cdfon response
figure(72);
cdfplot(allControlOn); hold on;
cdfplot(allL200On); hold on;
cdfplot(allL100On); hold on;
cdfplot(allL50On); hold on
cdfplot(allL10On); hold on;
cdfplot(allL1On); legend("Control","L200","L100","L50","L10","L1")
title('Cumulative Distribution all conditions On')
xlabel('time(ms)');ylabel('probability')
saveas(figure(72),[pwd '/New figures Jesus/all/cdfOnAll.pdf']);
saveas(figure(72),[pwd '/New figures Jesus/all/cdfOnAll.emf']);
%boxplots
allOff=[allControlOff,allL200Off,allL100Off,allL50Off,allL10Off,allL1Off];
allOfflabels=["Control","L200","L100","L50","L10","L1"];
figure(73); boxplot(allOff,allOfflabels,'Whisker',2)
title('Off response')
saveas(figure(73),[pwd '/New figures Jesus/all/Boxplot_AllOff.pdf']);
saveas(figure(73),[pwd '/New figures Jesus/all/Boxplot_AllOff.emf']);
allOn=[allControlOn,allL200On,allL100On,allL50On,allL10On,allL1On];
allOnlabels=["Control","L200","L100","L50","L10","L1"];
figure(74); boxplot(allOn,allOnlabels,'Whisker',2)
title('On Response')
saveas(figure(74),[pwd '/New figures Jesus/all/Boxplot_AllOn.pdf']);
saveas(figure(74),[pwd '/New figures Jesus/all/Boxplot_AllOn.emf']);


