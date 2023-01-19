figure(1);
subplot(3,6,[1 2 3])

MeanM1axonmap=mean(M1AxonMAP,1)
[sortS1axonmapNew idsort]=sort(MeanM1axonmap,'descend')
contour(M1AxonMAP(:,idsortA),'Fill','on','LevelStep',1)
view([90 90])
colormap(jet(128))
set(gca,'yTick',1:size(M1AxonMAP,1))
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:size(M1AxonMAP,1))
set(gca,'xTickLabel',rCaxonsMapNames(idsort))
set(gca,'yTick',1:13)
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:15)
set(gca,'xTickLabel',rCaxonsMapNames(idsortA))
title("MCRC rostro-caudal distribution axon map per structure")
n=max(M1AxonMAP(:))-1
n=n/5
colorbar('Ticks',[0 1 1+n 1+2*n 1+3*n 1+4*n 1+5*n ],'TickLabels',{'-','+','+','++','+++','++++','+++++'})
xlabel("Thalamic Nuclei")
hold on


subplot(3,6,[7 8 9 ])
MeanS1axonmap=mean(S1AxonMAPNew,1)
[sortS1axonmapNew idsort]=sort(MeanS1axonmap,'descend')
contour(S1AxonMAPNew(:,idsortA),'Fill','on','LevelStep',1)
view([90 90])
colormap(jet(128))
set(gca,'yTick',1:size(S1AxonMAPNew,1))
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:size(S1AxonMAPNew,1))
set(gca,'xTickLabel',rCaxonsMapNames(idsort))
set(gca,'yTick',1:13)
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:15)
set(gca,'xTickLabel',rCaxonsMapNames(idsortA))
title("BCRC rostro-caudal distribution axon map per structure")
n=max(S1AxonMAPNew(:))-1
n=n/5
colorbar('Ticks',[0 1 1+n 1+2*n 1+3*n 1+4*n 1+5*n ],'TickLabels',{'-','+','+','++','+++','++++','+++++'})
xlabel("Thalamic Nuclei")

hold on


subplot(3,6,[13 14 15])
MeanBSaxonmap=mean(BSAxonMAP,1)
[sortS1axonmapNew idsort]=sort(MeanBSaxonmap,'descend')
contour(BSAxonMAP(:,idsortA),'Fill','on','LevelStep',1)
view([90 90])
colormap(jet(128))
set(gca,'yTick',1:size(BSAxonMAP,1))
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:size(BSAxonMAP,1))
set(gca,'xTickLabel',rCaxonsMapNames(idsort))
set(gca,'yTick',1:13)
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:15)
set(gca,'xTickLabel',rCaxonsMapNames(idsortA))
title("BsRC rostro-caudal distribution axon map per structure")
n=max(BSAxonMAP(:))-1
n=n/5
colorbar('Ticks',[0 1 1+n 1+2*n 1+3*n 1+4*n 1+5*n ],'TickLabels',{'-','+','+','++','+++','++++','+++++'})
xlabel("Thalamic Nuclei")
hold off

%zscoreM1=zscore(mean(M1AxonMAP,1),0);
%[sortM1axonmapNew idsortA]=sort(zscoreM1,'descend');
%zscoreS1=zscore(mean(S1AxonMAPNew,1),0);
%zscoreBS=zscore(mean(BSAxonMAP,1),0);

orderMaxM1=mean(M1AxonMAP,1);
NMaxM1=orderMaxM1*1./max(orderMaxM1)
[SNMaxM1 idsortA]=sort(NMaxM1,'descend')

orderMaxS1=mean(S1AxonMAPNew,1);
NMaxS1=orderMaxS1*1./max(orderMaxS1)

orderMaxBS=mean(BSAxonMAP,1);
NMaxBS=orderMaxBS*1./max(orderMaxBS)


subplot(3,6,[5 6 11 12 16 17])

%plot(sortM1axonmapNew,'r--O','MarkerFaceColor','r')
%plot(SNMaxM1,'r--O','MarkerFaceColor','r')
%view([90 90])
%set(gca,'xTick',1:15)
%set(gca,'xTickLabel',rCaxonsMapNames(idsortA))
%hold on



%subplot(3,6,[5 6 11 12 16 17])
%plot(zscoreS1(idsortA),'b--O','MarkerFaceColor','b')
%plot(NMaxS1(idsortA),'b--O','MarkerFaceColor','b')
%view([90 90])
%set(gca,'xTick',1:15)
%set(gca,'xTickLabel',rCaxonsMapNames(idsortA))
hold on 


%zscoreALL=cat(2,zscoreM1(idsortA)',zscoreS1(idsortA)',zscoreBS(idsortA)')

NMaxALL=cat(2,SNMaxM1',NMaxS1(idsortA)',NMaxBS(idsortA)')
subplot(3,6,[5 6 11 12 16 17])

%plot(zscoreBS(idsortA),'g--O','MarkerFaceColor','g')
%plot(NMaxBS(idsortA),'g--O','MarkerFaceColor','g')
barall=bar(NMaxALL)
barall(1).FaceColor=[1 0 0];
hold on
barall(2).FaceColor=[0 0 1]
barall(3).FaceColor=[0 1 0]
view([90 90])
set(gca,'xTick',1:15)
set(gca,'xTickLabel',rCaxonsMapNames(idsortA)) 
set(gca,'yTick',0:0.25:1)
%set(gca,'yTickLabel',[0:0.25:1])
box off
%ylabel("zscore")
ylabel("Normalized to max. value")
ylim([0 1.5])
title("Comparison MC BC Bs")
legend("M1","S1","BS",'Location','southeast')