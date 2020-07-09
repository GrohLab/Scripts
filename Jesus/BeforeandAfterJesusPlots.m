%% % load files 
close all
%% get spike times L200 and plots condition L200
%do connected categorical scatter plots % add titles and labels
% add variables to change the desired window % measure median for each unit %add the box plot %add experimental id
%think about statistical test to see effect on single neuron level.
%next time with rebecca add dkl %need to find the appropriate nonparametric 2 factor anova
RspControlOn=ClustersControlOnset*1000; %now in ms 
RspControlOff=ClustersControlOnset*1000;
%find all unwanted spikes
excludefromOnset=[find(RspControlOn<2);find(RspControlOn>30)];
excludefromOffset=[find(RspControlOff<102);find(RspControlOff>130)];
RspControlOn(excludefromOnset)=nan;
RspControlOff(excludefromOffset)=nan;
yControlOn=nanmedian(RspControlOn);
yControlOff=nanmedian(RspControlOff);
%repeat for getting L200 
RspL200On=ClustersL200Onset*1000; % now in ms
RspL200Off=ClustersL200Onset*1000;
%find all unwanted spikes
excludefromOnsetL200=[find(RspL200On<2);find(RspL200On>30)];
excludefromOffsetL200=[find(RspL200Off<102);find(RspL200Off>130)];
RspL200On(excludefromOnsetL200)=nan;
RspL200Off(excludefromOffsetL200)=nan
yL200On=nanmedian(RspL200On);
yL200Off=nanmedian(RspL200Off);
x=ones(size(yL200On));
x1=ones(size(yL200On))*2;
% Population boxplot
%Off response
figure (1); subplot(1,2,1)
OffControlandL200=[yControlOff' yL200Off'];
OffLabels=["ControlOff","L200Off"]
boxplot(OffControlandL200,OffLabels);
title(["Off Rsp p= ",signrank(yControlOff,yL200Off)]);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL200Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2) %trends line per clusters
 OnControlandL200=[yControlOn' yL200On']
 OnLabels=["ControlOn","L200On"]
boxplot(OnControlandL200,OnLabels);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL200On(i)],'-o') 
end
xlim([0.75 2.25])
title(["On Rsp p=",signrank(yControlOn,yL200On)])
xlabel('Condition'); ylabel('time(ms)');
saveas(figure(1),[pwd '/new figures Jesus/ControlandL200/PoPboxplot On and Off Control and L200.emf']);
saveas(figure(1),[pwd '/new figures Jesus/ControlandL200/PoPboxplot On and Off Control and L200.pdf']);
%Boxplot individual clusters
%On respose
figure (2); subplot(1,2,1);
boxplot(RspControlOn,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control On Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL200On,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L200 On Response'),xlabel('times(ms)');ylabel('Clusters');
saveas(figure(2),[pwd '/New figures Jesus/ControlandL200/Boxplot individiual clasters On Control and L200.emf']);
saveas(figure(2),[pwd '/New figures Jesus/ControlandL200/Boxplot individiual clasters On Control and L200.pdf']);
%Off response
figure (3); subplot(1,2,1);
boxplot(RspControlOff,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control Off Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL200Off,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L200 Off Response');xlabel('times(ms)');ylabel('Clusters');
saveas(figure(3),[pwd '/New figures Jesus/ControlandL200/Boxplot individiual clasters Off Control and L200.emf']);
saveas(figure(3),[pwd '/New figures Jesus/ControlandL200/Boxplot individiual clasters Off Control and L200.pdf']);
% get all spikes On and Off for  scatterplots and cdf
allControlOn=RspControlOn(:); allL200On=RspL200On(:);
allControlOff=RspControlOff(:);allL200Off=RspL200Off(:);
% On scatter plot and Cumulative distribution
figure (4); scatterhist(allControlOn,allL200On,'Direction','out','NBins',50);
saveas(figure(4),[pwd '/New figures Jesus/ControlandL200/scatterhist On response Control and L200.emf']);
saveas(figure(4),[pwd '/New figures Jesus/ControlandL200/scatterhist On response Control and L200.pdf']);
figure (5); cdfplot(allControlOn)% off response scatter plot and Cumulative distribution
hold on
cdfplot(allL200On); title('Cumulative Distribution On Response'); legend('Control','L200','Location','southeast')
xlabel('time');ylabel('cumulative probability')
[plogicalOn,ksprobl200On]=kstest2(allControlOn,allL200On,'Alpha',0.05)
title(["Cumulative Distribution On Response p=",ksprobl200On]);
saveas(figure(5),[pwd '/New figures Jesus/ControlandL200/Cumulative Distribution On response Control and L200.emf']);
saveas(figure(5),[pwd '/New figures Jesus/ControlandL200/Cumulative Distribution On response Control and L200.pdf']);
%Off response scatter plot and cumulative distribution
figure (6); scatterhist(allControlOff,allL200Off,'Direction','out','NBins',50);
saveas(figure(6),[pwd '/New figures Jesus/ControlandL200/scatter plot off response control and L200 off.emf']);
saveas(figure(6),[pwd '/New figures Jesus/ControlandL200/scatter plot off response control and L200 off.pdf']);
figure (7); cdfplot(allControlOff); hold on;
cdfplot(allL200Off);  legend('Control','L200','Location','southeast');xlabel('time');ylabel('cumulative probability')
[plogicalOff,ksprobl200Off]=kstest2(allControlOff,allL200Off,'Alpha',0.05);
title(["Cumulative Distribution Off Response p=",ksprobl200Off]);
saveas(figure(7),[pwd '/New figures Jesus/ControlandL200/Cumulative distribution off response control and L200 off.emf']);
saveas(figure(7),[pwd '/New figures Jesus/ControlandL200/Cumulative distribution off response control and L200 off.pdf']);
%% %condition Laser 100 get spike times L100 and plots condition L100
%repeat for getting L100 
RspL100On=ClustersL100Onset*1000; % now in ms
RspL100Off=ClustersL100Onset*1000;
%find all unwanted spikes
excludefromOnsetL100=[find(RspL100On<2);find(RspL100On>30)];
excludefromOffsetL100=[find(RspL100Off<102);find(RspL100Off>130)];
RspL100On(excludefromOnsetL100)=nan;
RspL100Off(excludefromOffsetL100)=nan;
yL100On=nanmedian(RspL100On);
yL100Off=nanmedian(RspL100Off);
x=ones(size(yL100On));
x1=ones(size(yL100On))*2;
% Population boxplot control and L100
%Off response
figure (8); subplot(1,2,1);
OffControlandL100=[yControlOff' yL100Off'];
OffLabels100=["ControlOff","L100Off"]
boxplot(OffControlandL100,OffLabels100);
title(["Off Rsp ",signrank(yControlOff,yL100Off)]);
hold on % 
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL100Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2); %trends line per clusters
 OnControlandL100=[yControlOn' yL100On'];
 OnLabels100=["ControlOn","L100On"];
boxplot(OnControlandL100,OnLabels100);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL100On(i)],'-o') 
end
xlim([0.75 2.25])
title(["On Rsp",signrank(yControlOn,yL100On)])
xlabel('Condition'); ylabel('time(ms)');
saveas(figure(8),[pwd '/New figures Jesus/ControlandL100/PoPboxplot On and Off Control and L100.emf']);
saveas(figure(8),[pwd '/New figures Jesus/ControlandL100/PoPboxplot On and Off Control and L100.pdf']);
%Boxplot individual clusters
%On respose
figure (9); subplot(1,2,1)
boxplot(RspControlOn,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control On Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL100On,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L100 On Response'),xlabel('times(ms)');ylabel('Clusters');
saveas(figure(9),[pwd '/New figures Jesus/ControlandL100/Boxplot individiual clasters On Control and L100.emf']);
saveas(figure(9),[pwd '/New figures Jesus/ControlandL100/Boxplot individiual clasters On Control and L100.pdf']);
%Off response
figure (10); subplot(1,2,1);
boxplot(RspControlOff,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control Off Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL100Off,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L100 Off Response');xlabel('times(ms)');ylabel('Clusters');
saveas(figure(10),[pwd '/New figures Jesus/ControlandL100/Boxplot individiual clasters Off Control and L100.emf']);
saveas(figure(10),[pwd '/New figures Jesus/ControlandL100/Boxplot individiual clasters Off Control and L100.pdf']);
% get all spikes On and Off for  scatterplots and cdf
allControlOn=RspControlOn(:); allL100On=RspL100On(:);
allControlOff=RspControlOff(:);allL100Off=RspL100Off(:);
% On scatter plot and Cumulative distribution
figure (11); scatterhist(allControlOn,allL100On,'Direction','out','NBins',50);
saveas(figure(11),[pwd '/New figures Jesus/ControlandL100/scatterhist On response Control and L100.emf']);
saveas(figure(11),[pwd '/New figures Jesus/ControlandL100/scatterhist On response Control and L100.pdf']);
figure (12); cdfplot(allControlOn)% off response scatter plot and Cumulative distribution
hold on;
cdfplot(allL100On); title('Cumulative Distribution On Response'); legend('Control','L100','Location','southeast');
xlabel('time');ylabel('cumulative probability')
[plogicalOn100,ksprobl100On]=kstest2(allControlOn,allL100On,'Alpha',0.05);
title(["Cumulative Distribution On Response p=",ksprobl100On]);
saveas(figure(12),[pwd '/New figures Jesus/ControlandL100/Cumulative Distribution On response Control and L100.emf']);
saveas(figure(12),[pwd '/New figures Jesus/ControlandL100/Cumulative Distribution On response Control and L100.pdf']);
%Off response scatter plot and cumulative distribution
figure (13); scatterhist(allControlOff,allL100Off,'Direction','out','NBins',50);
saveas(figure(13),[pwd '/New figures Jesus/ControlandL100/scatter plot off response control and L100 off.emf']);
saveas(figure(13),[pwd '/New figures Jesus/ControlandL100/scatter plot off response control and L100 off.pdf']);
figure (14); cdfplot(allControlOff); hold on;
cdfplot(allL100Off);  legend('Control','L100','Location','southeast');xlabel('time');ylabel('cumulative probability')
[plogicalOff,ksprobl100Off]= kstest2(allControlOff,allL100Off,'Alpha',0.05);
title(["Cumulative Distribution Off Response p=",ksprobl100Off]);
saveas(figure(14),[pwd '/New figures Jesus/ControlandL100/Cumulative distribution off response control and L100 off.emf']);
saveas(figure(14),[pwd '/New figures Jesus/ControlandL100/Cumulative distribution off response control and L100 off.pdf']);
%% %Condition L50 get spike times and plots
%repeat for getting L50 
RspL50On=ClustersL50Onset*1000; % now in ms
RspL50Off=ClustersL100Onset*1000;
%find all unwanted spikes
excludefromOnsetL50=[find(RspL50On<2);find(RspL50On>30)];
excludefromOffsetL50=[find(RspL50Off<102);find(RspL50Off>130)];
RspL50On(excludefromOnsetL50)=nan;
RspL50Off(excludefromOffsetL50)=nan;
yL50On=nanmedian(RspL50On);
yL50Off=nanmedian(RspL50Off);
x=ones(size(yL50On));
x1=ones(size(yL50On))*2;
% Population boxplot control and L50
%Off response
figure (15); subplot(1,2,1);
OffControlandL50=[yControlOff' yL50Off'];
OffLabels50=["ControlOff","L50Off"];
boxplot(OffControlandL50,OffLabels50);
title(["Off Rsp ",signrank(yControlOff,yL50Off)]);
hold on % 
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL50Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2); %trends line per clusters
 OnControlandL50=[yControlOn' yL50On'];
 OnLabels50=["ControlOn","L50On"];
boxplot(OnControlandL50,OnLabels50);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL50On(i)],'-o');
end
xlim([0.75 2.25]);
title(["On Rsp",signrank(yControlOn,yL50On)]);
xlabel('Condition'); ylabel('time(ms)');
saveas(figure(15),[pwd '/New figures Jesus/ControlandL50/PoP boxplot On and Off Control and L50.emf']);
saveas(figure(15),[pwd '/New figures Jesus/ControlandL50/PoP boxplot On and Off Control and L50.pdf']);
%Boxplot individual clusters
%On respose
figure (16); subplot(1,2,1);
boxplot(RspControlOn,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control On Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL50On,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L50 On Response'),xlabel('times(ms)');ylabel('Clusters');
saveas(figure(16),[pwd '/New figures Jesus/ControlandL50/Boxplot individiual clasters On Control and L50.emf']);
saveas(figure(16),[pwd '/New figures Jesus/ControlandL50/Boxplot individiual clasters On Control and L50.pdf']);
%Off response
figure (17); subplot(1,2,1);
boxplot(RspControlOff,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control Off Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL50Off,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L50 Off Response');xlabel('times(ms)');ylabel('Clusters');
saveas(figure(17),[pwd '/New figures Jesus/ControlandL50/Boxplot individiual clasters Off Control and L50.emf']);
saveas(figure(17),[pwd '/New figures Jesus/ControlandL50/Boxplot individiual clasters Off Control and L50.pdf']);
% get all spikes On and Off for  scatterplots and cdf
allControlOn=RspControlOn(:); allL50On=RspL50On(:);
allControlOff=RspControlOff(:);allL50Off=RspL50Off(:);
% On scatter plot and Cumulative distribution
figure (18); scatterhist(allControlOn,allL50On,'Direction','out','NBins',50);
saveas(figure(18),[pwd '/New figures Jesus/ControlandL50/scatterhist On response Control and L50.emf']);
saveas(figure(18),[pwd '/New figures Jesus/ControlandL50/scatterhist On response Control and L50.pdf']);
figure (19); cdfplot(allControlOn)% off response scatter plot and Cumulative distribution
hold on
cdfplot(allL50On); title('Cumulative Distribution On Response'); legend('Control','L50','Location','southeast')
xlabel('time');ylabel('cumulative probability')
[plogicalOn50,ksprobl50On]=kstest2(allControlOn,allL50On,'Alpha',0.05);
title(["Cumulative Distribution On Response p=",ksprobl50On]);
saveas(figure(19),[pwd '/New figures Jesus/ControlandL50/Cumulative Distribution On response Control and L50.emf'])
saveas(figure(19),[pwd '/New figures Jesus/ControlandL50/Cumulative Distribution On response Control and L50.pdf'])%Off response scatter plot and cumulative distribution
figure (20); scatterhist(allControlOff,allL50Off,'Direction','out','NBins',50);
saveas(figure(20),[pwd '/New figures Jesus/ControlandL50/scatter plot off response control and L50 off.emf']);
saveas(figure(20),[pwd '/New figures Jesus/ControlandL50/scatter plot off response control and L50 off.pdf']);
figure (21); cdfplot(allControlOff); hold on;
cdfplot(allL50Off);  legend('Control','L50','Location','southeast');xlabel('time');ylabel('cumulative probability')
[plogicalOff,ksprobl50Off]= kstest2(allControlOff,allL50Off,'Alpha',0.05);
title(["Cumulative Distribution Off Response p=",ksprobl50Off]);
saveas(figure(21),[pwd '/New figures Jesus/ControlandL50/Cumulative distribution off response control and L50 off.emf']);
saveas(figure(21),[pwd '/New figures Jesus/ControlandL50/Cumulative distribution off response control and L50 off.pdf']);

%% %get spike times and plots Condition L10
%repeat for getting L10 
RspL10On=ClustersL10Onset*1000; % now in ms
RspL10Off=ClustersL10Onset*1000;
%find all unwanted spikes
excludefromOnsetL10=[find(RspL10On<2);find(RspL10On>30)];
excludefromOffsetL10=[find(RspL10Off<102);find(RspL10Off>130)];
RspL10On(excludefromOnsetL10)=nan;
RspL10Off(excludefromOffsetL10)=nan;
yL10On=nanmedian(RspL10On);
yL10Off=nanmedian(RspL10Off);
x=ones(size(yL10On));
x1=ones(size(yL10On))*2;
% Population boxplot control and L10
%Off response
figure (22); subplot(1,2,1);
OffControlandL10=[yControlOff' yL10Off'];
OffLabels10=["ControlOff","L10Off"];
boxplot(OffControlandL10,OffLabels10);
title(["Off Rsp ",signrank(yControlOff,yL10Off)]);
hold on % 
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL10Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2); %trends line per clusters
 OnControlandL10=[yControlOn' yL10On'];
 OnLabels10=["ControlOn","L10On"];
boxplot(OnControlandL10,OnLabels10);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL10On(i)],'-o');
end
xlim([0.75 2.25]);
title(["On Rsp",signrank(yControlOn,yL10On)]);
xlabel('Condition'); ylabel('time(ms)');
saveas(figure(22),[pwd '/New figures Jesus/ControlandL10/PoP boxplot On and Off Control and L10.emf']);
saveas(figure(22),[pwd '/New figures Jesus/ControlandL10/PoP boxplot On and Off Control and L10.pdf']);
%Boxplot individual clusters
%On respose
figure (23); subplot(1,2,1);
boxplot(RspControlOn,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control On Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL10On,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L10 On Response'),xlabel('times(ms)');ylabel('Clusters');
saveas(figure(23),[pwd '/New figures Jesus/ControlandL10/Boxplot individiual clasters On Control and L50.emf']);
saveas(figure(23),[pwd '/New figures Jesus/ControlandL10/Boxplot individiual clasters On Control and L50.pdf']);
%Off response
figure (24); subplot(1,2,1);
boxplot(RspControlOff,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control Off Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL10Off,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L10 Off Response');xlabel('times(ms)');ylabel('Clusters');
saveas(figure(24),[pwd '/New figures Jesus/ControlandL10/Boxplot individiual clasters Off Control and L10.emf']);
saveas(figure(24),[pwd '/New figures Jesus/ControlandL10/Boxplot individiual clasters Off Control and L10.pdf']);
% get all spikes On and Off for  scatterplots and cdf
allControlOn=RspControlOn(:); allL10On=RspL10On(:);
allControlOff=RspControlOff(:);allL10Off=RspL10Off(:);
% On scatter plot and Cumulative distribution
figure (25); scatterhist(allControlOn,allL10On,'Direction','out','NBins',50);
saveas(figure(25),[pwd '/New figures Jesus/ControlandL10/scatterhist On response Control and L10.emf']);
saveas(figure(25),[pwd '/New figures Jesus/ControlandL10/scatterhist On response Control and L10.pdf']);
figure (26); cdfplot(allControlOn)% off response scatter plot and Cumulative distribution
hold on
cdfplot(allL10On); title('Cumulative Distribution On Response'); legend('Control','L10','Location','southeast')
xlabel('time');ylabel('cumulative probability')
[plogicalOn10,ksprobl10On]=kstest2(allControlOn,allL10On,'Alpha',0.05);
title(["Cumulative Distribution On Response p=",ksprobl10On]);
saveas(figure(26),[pwd '/New figures Jesus/ControlandL10/Cumulative Distribution On response Control and L10.emf']);
saveas(figure(26),[pwd '/New figures Jesus/ControlandL10/Cumulative Distribution On response Control and L10.pdf']);
%Off response scatter plot and cumulative distribution
figure (27); scatterhist(allControlOff,allL10Off,'Direction','out','NBins',50);
saveas(figure(27),[pwd '/New figures Jesus/ControlandL10/scatter plot off response control and L10 off.emf']);
saveas(figure(27),[pwd '/New figures Jesus/ControlandL10/scatter plot off response control and L10 off.pdf']);
figure (28); cdfplot(allControlOff); hold on;
cdfplot(allL10Off);  legend('Control','L10','Location','southeast');xlabel('time');ylabel('cumulative probability')
[plogicalOff,ksprobl10Off]= kstest2(allControlOff,allL10Off,'Alpha',0.05);
title(["Cumulative Distribution Off Response p=",ksprobl10Off]);
saveas(figure(28),[pwd '/New figures Jesus/ControlandL10/Cumulative distribution off response control and L10 off.emf']);
saveas(figure(28),[pwd '/New figures Jesus/ControlandL10/Cumulative distribution off response control and L10 off.pdf']);

%% % get spikes times and plots condition L1

%repeat for getting L1
RspL1On=ClustersL1Onset*1000; % now in ms
RspL1Off=ClustersL1Onset*1000;
%find all unwanted spikes
excludefromOnsetL1=[find(RspL1On<2);find(RspL1On>30)];
excludefromOffsetL1=[find(RspL1Off<102);find(RspL1Off>130)];
RspL1On(excludefromOnsetL1)=nan;
RspL1Off(excludefromOffsetL1)=nan;
yL1On=nanmedian(RspL1On);
yL1Off=nanmedian(RspL1Off);
x=ones(size(yL1On));
x1=ones(size(yL1On))*2;
% Population boxplot control and L1
%Off response
figure (29); subplot(1,2,1);
OffControlandL1=[yControlOff' yL1Off'];
OffLabels1=["ControlOff","L1Off"];
boxplot(OffControlandL1,OffLabels1);
title(["Off Rsp ",signrank(yControlOff,yL1Off)]);
hold on % 
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL1Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2); %trends line per clusters
 OnControlandL10=[yControlOn' yL1On'];
 OnLabels1=["ControlOn","L1On"];
boxplot(OnControlandL10,OnLabels1);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL1On(i)],'-o');
end
xlim([0.75 2.25]);
title(["On Rsp",signrank(yControlOn,yL1On)]);
xlabel('Condition'); ylabel('time(ms)');
saveas(figure(29),[pwd '/New figures Jesus/ControlandL1/PoP boxplot On and Off Control and L1.emf']);
saveas(figure(29),[pwd '/New figures Jesus/ControlandL1/PoP boxplot On and Off Control and L1.pdf']);
%Boxplot individual clusters
%On respose
figure (30); subplot(1,2,1);
boxplot(RspControlOn,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control On Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL1On,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L1 On Response'),xlabel('times(ms)');ylabel('Clusters');
saveas(figure(30),[pwd '/New figures Jesus/ControlandL1/Boxplot individiual clasters On Control and L1.emf']);
saveas(figure(30),[pwd '/New figures Jesus/ControlandL1/Boxplot individiual clasters On Control and L1.pdf']);
%Off response
figure (31); subplot(1,2,1);
boxplot(RspControlOff,'Orientation','horizontal','BoxStyle','filled','Colors','k');
title('Clusters Control Off Response');
xlabel('times(ms)');ylabel('Clusters');
subplot(1,2,2);
boxplot(RspL1Off,'Orientation','horizontal','BoxStyle','filled','Colors','b');
title('Clusters L1 Off Response');xlabel('times(ms)');ylabel('Clusters');
saveas(figure(31),[pwd '/New figures Jesus/ControlandL1/Boxplot individiual clasters Off Control and L1.emf']);
saveas(figure(31),[pwd '/New figures Jesus/ControlandL1/Boxplot individiual clasters Off Control and L1.pdf']);
% get all spikes On and Off for  scatterplots and cdf
allControlOn=RspControlOn(:); allL1On=RspL1On(:);
allControlOff=RspControlOff(:);allL1Off=RspL1Off(:);
% On scatter plot and Cumulative distribution
figure (32); scatterhist(allControlOn,allL1On,'Direction','out','NBins',50);
saveas(figure(32),[pwd '/New figures Jesus/ControlandL1/scatterhist On response Control and L1.emf']);
saveas(figure(32),[pwd '/New figures Jesus/ControlandL1/scatterhist On response Control and L1.pdf']);
figure (33); cdfplot(allControlOn)% off response scatter plot and Cumulative distribution
hold on
cdfplot(allL10On); title('Cumulative Distribution On Response'); legend('Control','L1','Location','southeast');
xlabel('time');ylabel('cumulative probability');
[plogicalOn1,ksprobl1On]=kstest2(allControlOn,allL1On,'Alpha',0.05);
title(["Cumulative Distribution On Response p=",ksprobl1On]);
saveas(figure(33),[pwd '/New figures Jesus/ControlandL1/Cumulative Distribution On response Control and L1.emf']);
saveas(figure(33),[pwd '/New figures Jesus/ControlandL1/Cumulative Distribution On response Control and L1.pdf']);
%Off response scatter plot and cumulative distribution
figure (34); scatterhist(allControlOff,allL1Off,'Direction','out','NBins',50);
saveas(figure(34),[pwd '/New figures Jesus/ControlandL1/scatter plot off response control and L1 off.emf']),
saveas(figure(34),[pwd '/New figures Jesus/ControlandL1/scatter plot off response control and L1 off.pdf']),
figure (35); cdfplot(allControlOff); hold on;
cdfplot(allL10Off);  legend('Control','L1','Location','southeast');xlabel('time');ylabel('cumulative probability')
[plogicalOff,ksprobl1Off]= kstest2(allControlOff,allL1Off,'Alpha',0.05);
title(["Cumulative Distribution Off Response p=",ksprobl1Off]);
saveas(figure(35),[pwd '/New figures Jesus/ControlandL1/Cumulative distribution off response control and L1 off.emf']);
saveas(figure(35),[pwd '/New figures Jesus/ControlandL1/Cumulative distribution off response control and L1 off.pdf']);

