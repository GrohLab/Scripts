
%first subplot
figure(1);
ax1=subplot(13,4,[1 2 3 5 6 7 9 10 11 13 14 15 17 18 19 21 22 23 25 26 27 29 30 31])
contour(M1AxonMAP,'Fill','on','LevelStep',1)
view([90 90])
colormap(jet(128))
set(gca,'yTick',1:size(M1AxonMAP,1))
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:size(M1AxonMAP,1))
set(gca,'xTickLabel',rCaxonsMapNames)
set(gca,'yTick',1:13)
set(gca,'yTicklabel',[0.1:0.1:1.3])
set(gca,'xTick',1:17)
set(gca,'xTickLabel',rCaxonsMapNames)
title("M1RC rostro-caudal distribution axon map per structure")
n=max(M1AxonMAP(:))-1
n=n/5
colorbar('Ticks',[0 1 1+n 1+2*n 1+3*n 1+4*n 1+5*n ],'TickLabels',{'-','+','+','++','+++','++++','+++++'})
xlabel("Thalamic Nuclei")

%second subplot
ax2=subplot(13,4,[4 8 12 16 20 24 28 32])
bar(mean(M1AxonMAP,1))
view([90 90])
set(gca,'xTick',1:17)
set(gca,'xTickLabel',rCaxonsMapNames)
n1=max(mean(M1AxonMAP,1))
n1=n1/5
set(gca,'yTick',[0 1 1+n1 1+2*n1 1+3*n1 1+4*n1 1+5*n1])
set(gca,'yTicklabel',pluses)
title("M1RC Mean axon map per structure")
linkaxes([ax1, ax2],'x')

%third subplot
ax3=subplot(13,4,[37 38 39 41 42 43 45 46 47 49 50 51])
contour(M1AxonMAPBS,'Fill','on','LevelStep',1)
view([90 90])
set(gca,'xTick',1:7)
set(gca,'xTickLabel',rCaxonsMapBSNames)
set(gca,'yTick',1:size(M1AxonMAPBS,1))
set(gca,'yTicklabel',[0.1:0.1:(size(M1AxonMAPBS,1)/10)])
n3=max(M1AxonMAPBS(:))-1
n3=n3/5
colorbar('Ticks',[0 1 1+n3 1+2*n3 1+3*n3 1+4*n3 1+5*n3 ],'TickLabels',{'-','+','+','++','+++','++++','+++++'})
xlabel("BrainStem Nuclei")
ylabel("Rostro <-------> Caudal")
%4th figure
ax4=subplot(13,4,[40 44 48 52])
bar(mean(M1AxonMAPBS,1))
view([90 90])
set(gca,'xTick',1:7)
set(gca,'xTickLabel',rCaxonsMapBSNames)
n4=max(mean(M1AxonMAPBS,1))
n4=n4/5
set(gca,'yTick',[0 1 1+n4 1+2*n4 1+3*n4 1+4*n4 1+5*n4])
set(gca,'yTicklabel',pluses)
ylabel("Density")
linkaxes([ax3,ax4],'x')
saveas(figure(1),'\\lsdf02.urz.uni-heidelberg.de\sd19B001\Berin\Analysis\Matlab\RecipientApproach\AxonMaps\M1\All\M1ContourandMEanAllstructures.fig')




