figure
subplot(2,2,[1 3])

plot3((zeros(size(x3d(:,1),1),1)+0),flip(y3d(:,1))/0.61,flip(x3d(:,1))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
hold on

plot3((zeros(size(xIntS(:,1),1),1)+0),flip(yIntS)/0.61,flip(xIntS)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.8 0.8 0.8])
%plot3((zeros(size(xIntL(:,1),1),1)+0),flip(yIntL)/0.61,flip(xIntL)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(xIntI(:,1),1),1)+0),flip(yIntI)/0.61,flip(xIntI)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.8 0.8 0.8])

plot3((zeros(size(x3d(:,2),1),1)+50),flip(y3d(:,2))/0.61,flip(x3d(:,2))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,8),1),1)+100),flip(y3d(:,8))/0.61,flip(x3d(:,8))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,9),1),1)+150),flip(y3d(:,9))/0.61,flip(x3d(:,9))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,10),1),1)+200),flip(y3d(:,10))/0.61,flip(x3d(:,10))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,11),1),1)+250),flip(y3d(:,11))/0.61,flip(x3d(:,11))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,12),1),1)+300),flip(y3d(:,12))/0.61,flip(x3d(:,12))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,13),1),1)+350),flip(y3d(:,13))/0.61,flip(x3d(:,13))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,14),1),1)+400),flip(y3d(:,14))/0.61,flip(x3d(:,14))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,15),1),1)+450),flip(y3d(:,15))/0.61,flip(x3d(:,15))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,3),1),1)+500),flip(y3d(:,3))/0.61,flip(x3d(:,3))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,4),1),1)+550),flip(y3d(:,4))/0.61,flip(x3d(:,4))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,5),1),1)+600),flip(y3d(:,5))/0.61,flip(x3d(:,5))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,6),1),1)+650),flip(y3d(:,6))/0.61,flip(x3d(:,6))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,7),1),1)+700),flip(y3d(:,7))/0.61,flip(x3d(:,7))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1],'MarkerEdgeColor',[0.1 0.1 0.1])




box off
axis ij
xlabel("Antero-Posterior (µm)")
ylabel("Medio-lateral (µm)")
zlabel("dorso-Ventral (µm)")
title("3D distribution of Cortical and peripheral inputs in SC ")

pgreen1=plot3((zeros(size(x3d1(:,1),1),1)+0),flip(y3d1(:,1))/0.61,flip(x3d1(:,1))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
hold on
plot3((zeros(size(x3d1(:,2),1),2)+50),flip(y3d1(:,2))/0.61,flip(x3d1(:,2))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,8),1),1)+100),flip(y3d1(:,8))/0.61,flip(x3d1(:,8))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,9),1),1)+150),flip(y3d1(:,9))/0.61,flip(x3d1(:,9))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,10),1),1)+200),flip(y3d1(:,10))/0.61,flip(x3d1(:,10))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,11),1),1)+250),flip(y3d1(:,11))/0.61,flip(x3d1(:,11))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,12),1),1)+300),flip(y3d1(:,12))/0.61,flip(x3d1(:,12))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,13),1),1)+350),flip(y3d1(:,13))/0.61,flip(x3d1(:,13))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,14),1),1)+400),flip(y3d1(:,14))/0.61,flip(x3d1(:,14))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,15),1),1)+450),flip(y3d1(:,15))/0.61,flip(x3d1(:,15))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,3),1),1)+500),flip(y3d1(:,3))/0.61,flip(x3d1(:,3))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,4),1),2)+550),flip(y3d1(:,4))/0.61,flip(x3d1(:,4))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,5),1),1)+600),flip(y3d1(:,5))/0.61,flip(x3d1(:,5))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,6),1),1)+650),flip(y3d1(:,6))/0.61,flip(x3d1(:,6))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])
plot3((zeros(size(x3d1(:,7),1),1)+700),flip(y3d1(:,7))/0.61,flip(x3d1(:,7))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1],'MarkerEdgeColor',[0.1 1 0.1])

pred1=plot3((zeros(size(x3d2(:,1),1),1)+0),flip(y3d2(:,1))/0.61,flip(x3d2(:,1))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1])
hold on

plot3((zeros(size(x3d2(:,2),1),2)+50),flip(y3d2(:,2))/0.61,flip(x3d2(:,2))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,8),1),1)+100),flip(y3d2(:,8))/0.61,flip(x3d2(:,8))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,9),1),1)+150),flip(y3d2(:,9))/0.61,flip(x3d2(:,9))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,10),1),1)+200),flip(y3d2(:,10))/0.61,flip(x3d2(:,10))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,11),1),1)+250),flip(y3d2(:,11))/0.61,flip(x3d2(:,11))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,12),1),1)+300),flip(y3d2(:,12))/0.61,flip(x3d2(:,12))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,13),1),1)+350),flip(y3d2(:,13))/0.61,flip(x3d2(:,13))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,14),1),1)+400),flip(y3d2(:,14))/0.61,flip(x3d2(:,14))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,15),1),1)+450),flip(y3d2(:,15))/0.61,flip(x3d2(:,15))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,3),1),1)+500),flip(y3d2(:,3))/0.61,flip(x3d2(:,3))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,4),1),2)+550),flip(y3d2(:,4))/0.61,flip(x3d2(:,4))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,5),1),1)+600),flip(y3d2(:,5))/0.61,flip(x3d2(:,5))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,6),1),1)+650),flip(y3d2(:,6))/0.61,flip(x3d2(:,6))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])
plot3((zeros(size(x3d2(:,7),1),1)+700),flip(y3d2(:,7))/0.61,flip(x3d2(:,7))/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1],'MarkerEdgeColor',[1 0.1 0.1])

set(gca,'ZTick',0:1000:3000)
set(gca,'yTick',0:1000:3000)
legend([pgreen1 pred1],{'peripheral inputs','Cortical'},'Location','southeast')
[h icons]=legend([pgreen1 pred1],{'peripheral inputs','Cortical'},'Location','northeast')
icons(4).MarkerSize=5
icons(6).MarkerSize=5
text(15,500,1800,'Int','HorizontalAlignment','center')
text(15,500,2500,'Sup','HorizontalAlignment','center')
text(15,500,1100,'Deep','HorizontalAlignment','center')
view([-52 23]);
ylim([0.4 3000])
zticks(flip(mean(max(x3d,[],1))/0.61:-500:0))
set(gca,'ZTickLabel',flip([0:500:2712]))
yticks([0.1,1000,2000,3000])
set(gca,'YTickLabel',([0, 1000,2000,3000]))

ax1=subplot(2,2,2)
histogram(x3d1(:)./0.61,100,'FaceColor',[0 1 0],'Normalization','probability','EdgeColor',[0 1 0]) %pixel size =0.61
hold on
histogram(x3d2(:)./0.61,100,'FaceColor',[1 0.1 0.1],'Normalization','probability','EdgeColor',[1 0 0])
%set(gca,'XTick',0:500:2500)
xlabel("(µm)")
ylabel("probability")
title("peripheral inputs signal vs cortical input signal \newline                            ventral-dorsal")
xlabel("Ventral<---->Dorsal (µm)")
box off
le1=legend("Peripheral Inputs","Cortical Inputs",'Location','northwest')
le1.Box='off'
box off

xticks(flip(mean(max(x3d,[],1))/0.61:-500:0))
set(gca,'XTickLabel',flip([0:500:2712]))
set(gca,'XDir','reverse')
xlim(ax1, [0 mean(max(x3d,[],1))/0.61])

ax2=subplot(2,2,4)
histogram(y3d1(:)./0.61,100,'FaceColor',[0.1 1 0.1],'Normalization','probability','EdgeColor',[0 0.6 0])
hold on
histogram(y3d2(:)./0.61,100,'FaceColor',[1 0.1 0.1],'Normalization','probability','EdgeColor',[0.6 0 0])
set(gca,'XTick',0:500:2500)
xlabel("medio-lateral (µm)")
ylabel("probability")
%title("medio-lateral (µm)")
box off
le2=legend("Peripheral Inputs","Cortical Inputs",'Location','northwest')
le2.Box='off'

xlim(ax2, [0 mean(max(y3d,[],1))/0.61])
linkaxes([ax1 ax2],'x')
%xticks(flip(mean(max(y3d,[],1))/0.61:-500:0))


%% 
figure

ax0=subplot(4,4,[2 3 4 6 7 8 10 11 12])

plot(y3d./0.61,x3d./0.61,'k*','MarkerSize',0.5)
ylim([0 3000])
hold on
plot(y3d1(:)/0.61,x3d1(:)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 1 0.1])
plot(y3d2(:)/0.61,x3d2(:)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[1 0.1 0.1])
plot(flip(yIntS)/0.61,flip(xIntS)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.8 0.8 0.8])
%plot3((zeros(size(xIntL(:,1),1),1)+0),flip(yIntL)/0.61,flip(xIntL)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot(flip(yIntI)/0.61,flip(xIntI)/0.61,'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.8 0.8 0.8])

text(500,1800,'Int','HorizontalAlignment','center')
text(500,2500,'Sup','HorizontalAlignment','center')
text(500,1100,'Deep','HorizontalAlignment','center')
title("D-V and M-L distribution of Cortical and peripheral inputs in SC ")
yticks(flip(mean(max(x3d,[],1))/0.61:-500:0))
set(gca,'yTickLabel',flip([0:500:2712]))
% Obvious effect ;-)
% set(get(ax0, "XAxis"), "Visible", "off")
% set(get(ax0, "YAxis"), "Visible", "off")

ax1=subplot(4,4,[1 5 9])

histogram(x3d1(:)./0.61,100,'FaceColor',[0.1 1 0.1],'Normalization','probability','EdgeColor',[0.1 1 0.1],'Orientation', 'horizontal')
hold on
histogram(x3d2(:)./0.61,100,'FaceColor',[1 0.1 0.1],'Normalization','probability','EdgeColor',[1 0.1 0.1],'Orientation', 'horizontal')
set(gca,'YTick',0:500:2500)
ylabel("(µm)")
xlabel("probability")
title("ventral-dorsal")


box off
le1=legend("Peripheral Inputs","Cortical Inputs",'Location','northwest')
le1.Box='off'
box off
% view([-90 90])
ylim([0 3000])
yticks(flip(mean(max(x3d,[],1))/0.61:-500:0))
set(gca,'yTickLabel',flip([0:500:2712]))
linkaxes([ax0, ax1], 'y')
set(ax1, 'XDir',  'reverse')

ax2=subplot(4,4,[14 15 16])
histogram(y3d1(:)./0.61,100,'FaceColor',[0.1 1 0.1],'Normalization','probability','EdgeColor',[0.1 1 0.1])
hold on
histogram(y3d2(:)./0.61,100,'FaceColor',[1 0.1 0.1],'Normalization','probability','EdgeColor',[1 0.1 0.1])
set(gca,'XTick',0:500:2500)
xlabel("(µm)")
ylabel("probability")
title("medio-lateral (µm)")
box off
le2=legend("Peripheral Inputs","Cortical Inputs",'Location','northwest')
le2.Box='off'
xlim([0 3000])
% set(ax2, 'YDir',  'reverse')
set(gcf, 'Color', 'w')
set([ax0, ax1, ax2], 'Box','off','Color','none')
%%test
%hold on
[c1 c2]=histcounts((x3d1(:)./0.61),100)
c3=(c1./(sum(c1(:))))
c3=c3.*3000/max(c3)

[p1 p2]=histcounts((x3d2(:)./0.61),100)
p3=(p1./(sum(p1(:))))
p3=p3.*3000/max(p3)

plot(c3,c2(1:end-1))
hold on
plot(p3,p2(1:end-1))
ylim([0 3000])
%plot(c2(1:end-1),ans)

