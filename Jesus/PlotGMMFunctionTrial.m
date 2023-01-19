function PlotGMMSC(CordMatCVG,CordMatFDIO,XtestC,YtestC)

%[XtestC YtestC]=find(testGM2>=200);
figure
gmtestCVGMC=fitgmdist([CordMatCVG(:,1)/0.61,CordMatCVG(:,2)/0.61],2)
DM=fcontour(@(x,y) pdf(gmtestCVGMC, [y,x]), [[00, 00]/0.61;[2000,2000]],'LineWidth',1.5)
hold on

maxDenM=max(DM.ZData(:))
[XMaxTestDM YMaxTestDM]=find(DM.ZData==maxDenM)
NcellM=size(DM.ZData,1);rxValueM=DM.XRange(2)/NcellM
plot(CordMatFDIO(:,2)/0.61,CordMatFDIO(:,1)/0.61,'ro','MarkerSize',2,'MarkerEdgeColor',[0.7 0.2 0.2],'MarkerFaceColor',[0.7 0.2 0.2])
hold on
plot(XtestC/.61,YtestC/0.61,'ko','MarkerSize',1)
plot(CordMatCVG(:,2)/0.61,CordMatCVG(:,1)/0.61,'go','MarkerSize',2,'MarkerEdgeColor',[0.2 0.5 0.2],'MarkerFaceColor',[0.2 0.5 0.2])
plot(YMaxTestDM*rxValueM,XMaxTestDM*rxValueM,'Ko','MarkerSize',10,'MarkerFaceColor',[0.2 0.2 0.2])
xlim([2 2500])
ylim([2 2500])
view([90 90])
legend('GMM','RN','slice Contour','CVG-RN','Max prob CVG area','Location','southwest')
legend('boxoff')
title(" gausians mixture model MC CVG-RN")
ylabel("Dorso-Ventral")
xlabel("Dorso-Ventral")
ylabel("Medio-Lateral")
colormap(jet(128))
annotation('line',[0.618 0.618],[0.56 0.11],'LineWidth',1,'LineStyle','--');
annotation('line',[0.13 0.618],[0.585 0.585],'LineWidth',1,'LineStyle','--');
%pie
axes('Position',[0.625 0.625 .3 .3])
TestPieBC=pie([(length(CordMatFDIO)-length(CordMatCVG)) length(CordMatCVG)])
TestPieBC(1).FaceColor=[0.7 0.1 0.1]
TestPieBC(3).FaceColor=[0.1 0.7 0.1]
TestPieBC(4).HorizontalAlignment='left'
TestPieBC(2).HorizontalAlignment='right'
lgGMMMC=legend("RN","CVG-RN",'Location','south')
lgGMMMC.Position= [0.72    0.57    0.1639    0.0717]
legend('boxoff')