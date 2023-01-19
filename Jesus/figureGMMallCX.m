figure


%plot(XtestFB,YtestFB,'o','MarkerSize',1,'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3])


%plot(XtestCVB,YtestCVB,'go','MarkerSize',1)
gmtestCVGBC=fitgmdist([XYBCGMMCVGALLCX(:,1)/0.61,XYBCGMMCVGALLCX(:,2)/0.61],2)
D=fcontour(@(x,y) pdf(gmtestCVGALL, [y,x]), [[00, 00]/0.61;[2000,2000]],'LineWidth',1.5,'LevelStep',0.0000001,'Fill','on')
hold on

maxDen=max(D.ZData(:))
[XMaxTestD YMaxTestD]=find(D.ZData==maxDen)
Ncell=size(D.ZData,1);rxValue=D.XRange(2)/Ncell
plot(YMaxTestD*rxValue,XMaxTestD*rxValue,'Ko','MarkerSize',10,'MarkerFaceColor',[0.2 0.2 0.2])
plot(XYGMMFDIOALLCX(:,2)/0.61,XYGMMFDIOALLCX(:,1)/0.61,'ro','MarkerSize',2,'MarkerEdgeColor',[0.7 0.2 0.2],'MarkerFaceColor',[0.7 0.2 0.2])
hold on
plot(XtestC/.61,YtestC/0.61,'ko','MarkerSize',1)
plot(XYBCGMMCVGALLCX(:,2)/0.61,XYBCGMMCVGALLCX(:,1)/0.61,'go','MarkerSize',2,'MarkerEdgeColor',[0.2 0.5 0.2],'MarkerFaceColor',[0.2 0.5 0.2])
xlim([2 2500])
ylim([2 2500])
view([90 90])
legend('GMM','Max prob CVG area','RN','slice Contour','CVG-RN','Location','southwest')
legend('boxoff')
title(" gausians mixture model all cx CVG-RN area")
ylabel("Dorso-Ventral")
xlabel("Dorso-Ventral")
ylabel("Medio-Lateral")
colormap(jet(128))