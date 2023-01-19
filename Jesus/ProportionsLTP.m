figure

ColorPie=[0.8 0.8 0.8;0 0 0;0.4 0.4 0.4]

bar([0 4],[13 15],0.3,'FaceColor',[0 0 1])
hold on
bar([1 5],[21 10],0.3,'FaceColor',[1 0 0])
set(gca,'xTick',[0.5 4.5]);
set(gca,'xTickLabel',["Potentiated" "Depressed"])
set(gca,'YTick',[0 10 20])
ylabel("% wrt total responsive "); title ("Percentage of Modulated Clusters")
xlim([-3 14]);ylim([0 30])
text([10.5], [28],'L6','HorizontalAlignment','center','fontsize',13)
text([10.5], [13],'Sham','HorizontalAlignment','center','fontsize',13)
box off
annotation('line',[0.26 0.32],[0.73 0.73],'LineWidth',1);
text([0.5], [24],'p=0.011','HorizontalAlignment','center','fontsize',9)
annotation('line',[0.44 0.5],[0.60 0.60],'LineWidth',1);
text([4.5], [19],'p=0.12','HorizontalAlignment','center','fontsize',9)
legend('Sham Induction','L6 Induction')
set(legend,'Position',[0.3 0.85 0.05 0.05])
%insets1
axes('Position',[0.60 0.55 0.27 0.27])
labelsL6={'\leftarrow\rightarrow(68%) ','\uparrow(21%)','\downarrow(10%)'}
P1=pie([217 68 33],[0 1 1],labelsL6)
colormap(ColorPie)
P1(6).Position=[0.9 1.3 0]
P1(6).FontWeight='bold'
P1(4).Position=[1.5 0.7 0]
P1(4).FontWeight='bold'
P1(2).Position=[-0.0 -0.4 0]
hold off

box off 
%inset 2
axes('Position',[0.6 0.15 0.27 0.27])
labelsMock={'\leftarrow\rightarrow(72%) ','\uparrow(13%)','\downarrow(15%)'}
P2=pie([180 33 38],[0 1 1],labelsMock)

P2(6).Position=[0.85 1.2 0]
P2(6).FontWeight='bold'
P2(4).Position=[1.5 0.5 0]
P2(4).FontWeight='bold'
P2(2).Position=[-0.0 -0.4 0]

box off