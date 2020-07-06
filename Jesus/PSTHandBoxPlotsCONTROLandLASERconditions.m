
mkdir ('New figures Jesus/ControlandL200'); mkdir ('New figures Jesus/ControlandL100');
mkdir ('New figures Jesus/ControlandL50'); mkdir ('New figures Jesus/ControlandL10');
mkdir ('New figures Jesus/ControlandL1');

%% boxplot and psth control and L200
nbin = 350

figure(40);
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL200Onset,'Orientation','horizontal')
title('L200')
saveas(figure(40),[pwd '/New figures Jesus/ControlandL200/boxplotclustersOnandOffl200.emf']);
saveas(figure(40),[pwd '/New figures Jesus/ControlandL200/boxplotclustersOnandOffl200.pdf']);

figure (41) %psth for conditions LASER 200 y Control
allControl=ClustersControlOnset(:);
allL200=ClustersL200Onset(:);
subplot(3,1,1) %psth control
histogram(allControl,nbin)
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
strOn={'piezo on'}
text(0.1,390,strOff)
text(0,390,strOn)
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
subplot(3,1,2)% l200 condition
histogram(allL200,nbin,'FaceColor','r')
xlabel('time(s)'),ylabel('counts')
title('Population Psth L200')
hold on
plot([-0.200;-0.200], [0;400],'--','LineWidth',2,'Color','b')
hold on
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-0.2,350,strLaserOn)
text(0.1,370,strLaserOff)
text(0.1,390,strOff)
text(-0.0,390,strOn)
subplot(3,1,3)%overlapped histograms
histogram(allControl,350,'EdgeColor','b','DisplayStyle','stairs')
hold on
histogram(allL200,350,'DisplayStyle','stairs','EdgeColor','r','FaceAlpha',1)
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
text(0.1,390,strOff)
text(-0.0,390,strOn)
xlabel('time(s)'),ylabel('counts')
title('Population Psth L200')
legend ('Pop Psth Control', 'PoP Psth L200','Location','southwest')
saveas(figure(41),[pwd '/New figures Jesus/ControlandL200/Psth Control andL200.pdf']);
saveas(figure(41),[pwd '/New figures Jesus/ControlandL200/Psth Control andL200.emf']);
%% Boxplot and psth control and l100 condition
figure (42)
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL100Onset,'Orientation','horizontal')
title('L100')
saveas(figure(42),[pwd '/New figures Jesus/ControlandL100/boxplotclustersOnandOffl100.emf']);
saveas(figure(42),[pwd '/New figures Jesus/ControlandL100/boxplotclustersOnandOffl100.pdf']);

figure (43) %psth for conditions LASER 100 y Control
allControl=ClustersControlOnset(:);
allL100=ClustersL100Onset(:);
subplot(3,1,1) %psth control
histogram(allControl,nbin)
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
strOn={'piezo on'}
text(0.1,390,strOff)
text(0,390,strOn)
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
subplot(3,1,2)% l100 condition
histogram(allL100,nbin,'FaceColor','r')
xlabel('time(s)'),ylabel('counts')
title('Population Psth L100')
hold on
plot([-0.100;-0.100], [0;400],'--','LineWidth',2,'Color','b')
hold on
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-0.1,350,strLaserOn)
text(0.1,370,strLaserOff)
text(0.1,390,strOff)
text(-0.0,390,strOn)
subplot(3,1,3)%overlapped histograms
histogram(allControl,350,'EdgeColor','b','DisplayStyle','stairs')
hold on
histogram(allL100,350,'DisplayStyle','stairs','EdgeColor','r','FaceAlpha',1)
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
text(0.1,390,strOff)
text(-0.0,390,strOn)
xlabel('time(s)'),ylabel('counts')
title('Population Psth L100')
legend ('Pop Psth Control', 'PoP Psth L100','Location','southwest')
saveas(figure(43),[pwd '/New figures Jesus/ControlandL100/Psth Control andL100.pdf']);
saveas(figure(43),[pwd '/New figures Jesus/ControlandL100/Psth Control andL100.emf']);
%% boxplot and psth control and L50
figure (44)
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL50Onset,'Orientation','horizontal')
title('L50')
saveas(figure(44),[pwd '/New figures Jesus/ControlandL50/boxplotclustersOnandOffl50.emf']);
saveas(figure(44),[pwd '/New figures Jesus/ControlandL50/boxplotclustersOnandOffl50.pdf']);

figure (45) %psth for conditions LASER 50 y Control
allControl=ClustersControlOnset(:);
allL50=ClustersL50Onset(:);
subplot(3,1,1) %psth control
histogram(allControl,nbin)
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
strOn={'piezo on'}
text(0.1,390,strOff)
text(0,390,strOn)
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
subplot(3,1,2)% l100 condition
histogram(allL50,nbin,'FaceColor','r')
xlabel('time(s)'),ylabel('counts')
title('Population Psth L100')
hold on
plot([-0.05;-0.05], [0;400],'--','LineWidth',2,'Color','b')
hold on
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-0.1,350,strLaserOn)
text(0.1,370,strLaserOff)
text(0.1,390,strOff)
text(-0.0,390,strOn)
subplot(3,1,3)%overlapped histograms
histogram(allControl,350,'EdgeColor','b','DisplayStyle','stairs')
hold on
histogram(allL50,350,'DisplayStyle','stairs','EdgeColor','r','FaceAlpha',1)
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
text(0.1,390,strOff)
text(-0.0,390,strOn)
xlabel('time(s)'),ylabel('counts')
title('Population Psth L50')
legend ('Pop Psth Control', 'PoP Psth L50','Location','southwest')
saveas(figure(45),[pwd '/New figures Jesus/ControlandL50/Psth Control andL50.pdf']);
saveas(figure(45),[pwd '/New figures Jesus/ControlandL50/Psth Control andL50.emf']);
%% Boxplot and Psth Control and L10 condition

figure (46)
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL50Onset,'Orientation','horizontal')
title('L10')
saveas(figure(46),[pwd '/New figures Jesus/ControlandL10/boxplotclustersOnandOffl10.emf']);
saveas(figure(46),[pwd '/New figures Jesus/ControlandL10/boxplotclustersOnandOffl10.pdf']);

figure (47) %psth for conditions LASER 10 y Control

allControl=ClustersControlOnset(:);
allL10=ClustersL10Onset(:);
subplot(3,1,1) %psth control
histogram(allControl,nbin)
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
strOn={'piezo on'}
text(0.1,390,strOff)
text(0,390,strOn)
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
subplot(3,1,2)% l100 condition
histogram(allL10,nbin,'FaceColor','r')
xlabel('time(s)'),ylabel('counts')
title('Population Psth L100')
hold on
plot([-0.01;-0.01], [0;400],'--','LineWidth',2,'Color','b')
hold on
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-0.1,350,strLaserOn)
text(0.1,370,strLaserOff)
text(0.1,390,strOff)
text(-0.0,390,strOn)
subplot(3,1,3)%overlapped histograms
histogram(allControl,350,'EdgeColor','b','DisplayStyle','stairs')
hold on
histogram(allL10,350,'DisplayStyle','stairs','EdgeColor','r','FaceAlpha',1)
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
text(0.1,390,strOff)
text(-0.0,390,strOn)
xlabel('time(s)'),ylabel('counts')
title('Population Psth L10')
legend ('Pop Psth Control', 'PoP Psth L10','Location','southwest')
saveas(figure(47),[pwd '/New figures Jesus/ControlandL10/Psth Control andL10.pdf']);
saveas(figure(47),[pwd '/New figures Jesus/ControlandL10/Psth Control andL10.emf']);
%% Boxplot and Psth control and L1 condition

figure (48)
subplot(1,2,1)
boxplot(ClustersControlOnset,'Orientation','horizontal')
title('Control')
xlabel('time(s)'), ylabel('Clusters')
subplot(1,2,2)
boxplot(ClustersL1Onset,'Orientation','horizontal')
title('L1')
saveas(figure(48),[pwd '/New figures Jesus/ControlandL1/boxplotclustersOnandOffl1.emf']);
saveas(figure(48),[pwd '/New figures Jesus/ControlandL1/boxplotclustersOnandOffl1.pdf']);

figure (49) %psth for conditions LASER 1 y Control

allControl=ClustersControlOnset(:);
allL1=ClustersL1Onset(:);
subplot(3,1,1) %psth control
histogram(allControl,nbin)
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
strOn={'piezo on'}
text(0.1,390,strOff)
text(0,390,strOn)
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
subplot(3,1,2)% l100 condition
histogram(allL1,nbin,'FaceColor','r')
xlabel('time(s)'),ylabel('counts')
title('Population Psth L100')
hold on
plot([-0.001;-0.001], [0;400],'--','LineWidth',2,'Color','b')
hold on
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-0.1,350,strLaserOn)
text(0.1,370,strLaserOff)
text(0.1,390,strOff)
text(-0.0,390,strOn)
subplot(3,1,3)%overlapped histograms
histogram(allControl,350,'EdgeColor','b','DisplayStyle','stairs')
hold on
histogram(allL1,350,'DisplayStyle','stairs','EdgeColor','r','FaceAlpha',1)
plot([0.0;-0.0], [0;400],'--','LineWidth',2,'Color','g')
plot([0.1;0.1], [0;400],'--','LineWidth',2,'Color','g')
text(0.1,390,strOff)
text(-0.0,390,strOn)
xlabel('time(s)'),ylabel('counts')
title('Population Psth L1')
legend ('Pop Psth Control', 'PoP Psth L1','Location','southwest');
saveas(figure(49),[pwd '/New figures Jesus/ControlandL1/Psth Control andL1.pdf']);
saveas(figure(49),[pwd '/New figures Jesus/ControlandL1/Psth Control andL1.emf']);