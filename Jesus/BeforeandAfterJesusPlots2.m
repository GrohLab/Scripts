%% % load files 
%load('jittering_2mw_multiprobeE1-1_exportSpkTms.mat')
mkdir ('New figures Jesus/ControlandL200'); mkdir ('New figures Jesus/ControlandL100');
mkdir ('New figures Jesus/ControlandL50'); mkdir ('New figures Jesus/ControlandL10');
mkdir ('New figures Jesus/ControlandL1');mkdir('New figures Jesus/all')
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
yControlstdOn=nanstd(RspControlOn);
%repeat for getting L200 
yControlOff=nanmedian(RspControlOff);
yControlstdOff=nanstd(RspControlOff);
RspL200On=ClustersL200Onset*1000; % now in ms
RspL200Off=ClustersL200Onset*1000;
%find all unwanted spikes
excludefromOnsetL200=[find(RspL200On<2);find(RspL200On>30)];
excludefromOffsetL200=[find(RspL200Off<102);find(RspL200Off>130)];
RspL200On(excludefromOnsetL200)=nan;
RspL200Off(excludefromOffsetL200)=nan;
yL200On=nanmedian(RspL200On);
yL200stdOn=nanstd(RspL200On);
yL200Off=nanmedian(RspL200Off);
yL200stdOff=nanstd(RspL200Off);
x=ones(size(yL200On));
x1=ones(size(yL200On))*2;
% Population boxplot
%Off response
figure (1); subplot(1,2,1)
OffControlandL200=[yControlOff' yL200Off'];
OffControlandL200std=[yControlstdOff yL200stdOff];
OffLabels=["ControlOff","L200Off"]
boxplot(OffControlandL200,OffLabels);
title(["Median Off Rsp p= ",signrank(yControlOff,yL200Off)]);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOff(i) yL200Off(i)],'-o') ;
end
xlim([0.75 2.25]);
xlabel('Condition'); ylabel('time(ms)');
% on response
subplot(1,2,2) %trends line per clusters
 OnControlandL200=[yControlOn' yL200On']
 OnControlandL200std=[yControlstdOn' yL200stdOn'];
 OnLabels=["ControlOn","L200On"];
 boxplot(OnControlandL200,OnLabels);
hold on
for i=1:numel(x)
   plot([1.25 1.75],[yControlOn(i) yL200On(i)],'-o') 
end
xlim([0.75 2.25])
title(["Median On Rsp p=",signrank(yControlOn,yL200On)])
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

hold off;
iqrL200On=iqr(RspL200On);
iqrL200Off=iqr(RspL200Off);
iqrControlOn=iqr(RspControlOn);
iqrControlOff=iqr(RspControlOff);
figure;  boxplot(OnControlandL200,OnLabels);
figure;boxplot(OnControlandL200std,OnLabels);
title(["Std On Rsp p=",signrank(iqrL200Off,iqrControlOff)])


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
figure (12); scatterhist(allControlOff,allL100Off,'Direction','out','NBins',50);
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
figure (32); scatterhist(allControlOn,allL1On,'Direction','out','NBins',60);
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

%%
%Contour plots
%On response
XValidControlOn=allControlOn(~isnan(allControlOn));
YValidL1On=allL1On(~isnan(allL1On));
allValidOn=nan(4000,2);
allValidOn(1:size(XValidControlOn),2)=XValidControlOn;
allValidOn(1:size(YValidL1On),1)=YValidL1On;
HistAllOn=hist3(allValidOn,[30 30]);
figure(119);
contour(HistAllOn,'fill','on');
xlabel('spike times Control');ylabel(' spike times Laser 1');
colorbar; colormap(jet(40));
title('On Reponse');
figure(120); surf(hist3(allValidOn,[30 30]),'FaceAlpha',1);view([115,54])
colorbar; colormap(jet(10));
xlabel('Spike Times Control');ylabel('Spike Times Laser 200');
title('On Reponse');
colorbar; colormap(jet(40));
saveas(figure(120),[pwd '/New figures Jesus/ControlandL200/Surf ON ControlsVSlaser200.emf']);
saveas(figure(120),[pwd '/New figures Jesus/ControlandL200/Surf ON ControlsVSlaser200.pdf']);
saveas(figure(119),[pwd '/New figures Jesus/ControlandL200/Contour ON ControlsVSlaser200.emf']);
saveas(figure(119),[pwd '/New figures Jesus/ControlandL200/Contour ON ControlsVSlaser200.pdf']);

%% Contours
%Off response
XValidControlOff=allControlOff(~isnan(allControlOff));
YValidL1Off=allL1Off(~isnan(allL1Off));
allValidOff=nan(4000,2);
allValidOff(1:size(XValidControlOff),2)=XValidControlOff;
allValidOff(1:size(YValidL200Off),1)=YValidL200Off;
HistAllOff=hist3(allValidOff,[30 30]);
figure (122);
contour(HistAllOff,'fill','on');
colorbar; colormap(jet(40));
xlabel('Control');ylabel('Laser');
title('Off Reponse');
figure(123); surf(hist3(allValidOff,[30 30]),'FaceAlpha',1); view([115,54])
colorbar; colormap(jet(10));
xlabel('Spike times Control');ylabel('spike times Laser 1');
title('Off Reponse');
colorbar; colormap(jet(40));
saveas(figure(123),[pwd '/New figures Jesus/ControlandL200/Surf OFF ControlsVSlaser200.emf']);
saveas(figure(123),[pwd '/New figures Jesus/ControlandL200/Surf OFF ControlsVSlaser200.pdf']);
saveas(figure(122),[pwd '/New figures Jesus/ControlandL200/Contour OFF ControlsVSlaser200.emf']);
saveas(figure(122),[pwd '/New figures Jesus/ControlandL200/Contour OFF ControlsVSlaser200.pdf']);
%%
a=yL200On-yControlOn
BiggerMedianClusters=a>0.00;
SmallerMedianClusters=a<0.00;

ReducingMedianClustersControl=ClustersControlOnset(:,SmallerMedianClusters);
ReducingMedianClustersControl=ReducingMedianClustersControl;

IncreasingMedianClustersControl=ClustersControlOnset(:,BiggerMedianClusters);
AllReducingMedianControl=ReducingMedianClustersControl(:);
AllReducingMedianControlms=AllReducingMedianControl*1000;
figure; histogram(AllReducingMedianControlms,350,'EdgeColor','b','displayStyle','stairs','linewidth',1);
hold on
ReducingMedianClustersL200=ClustersL200Onset(:,SmallerMedianClusters);
ReducingMedianClustersL200=ReducingMedianClustersL200;
 
IncreasingMedianClustersL200=ClustersL200Onset(:,BiggerMedianClusters);
AllReducingMedianL200=ReducingMedianClustersL200(:);
AllReducingMedianL200ms=AllReducingMedianL200*1000;
histogram(AllReducingMedianL200ms,350,'EdgeColor','r','DisplayStyle','stairs','linewidth',1); hold on;
xlabel('time(s)'),ylabel('counts')
title('Population Psth Control')
hold on
plot([0.0;0.0], [0;400],'--','LineWidth',2,'Color','g')
strOn={'piezo on'}
text(0,370,strOn)
hold on;
plot([100;100],[0;400],'--','LineWidth',2,'Color','g')
strOff={'piezo off'}
text(100,370,strOff); hold on
plot([-200;-200], [0;400],'--','LineWidth',2,'Color','b')
strLaserOn={'Laser On'}
strLaserOff={'Laser Off'}
text(-200,370,strLaserOn)
text(100,360,strLaserOff)
legend ('Pop Psth Control', 'PoP Psth L200','Location','northwest')

%% boxplot reducing medians
figure(110)
subplot(1,2,1)
boxplot(RspControlOn(:,SmallerMedianClusters),'Orientation','horizontal','BoxStyle','filled')
subplot(1,2,2)
boxplot(RspL200On(:,SmallerMedianClusters),'Orientation','horizontal','Colors','r','BoxStyle','filled')

figure(111)
subplot(1,2,1)
boxplot(RspControlOff(:,SmallerMedianClusters),'Orientation','horizontal','BoxStyle','filled')
subplot(1,2,2)
boxplot(RspL200Off(:,SmallerMedianClusters),'Orientation','horizontal','Colors','r','BoxStyle','filled')

%% contour reducing medians ON
ReducingClustersRMControl=RspControlOn(:,SmallerMedianClusters);
allRMControlOn=ReducingClustersRMControl(:);
ReducingClustersRML200=RspL200On(:,SmallerMedianClusters);
allRML200On=ReducingClustersRML200(:);
boxplot([allRMControlOn allRML200On],OnLabels)
title(["On Response p=" ,signrank(allRMControlOn,allRML200On)])
XValidRMControlOn=allRMControlOn(~isnan(allRMControlOn));
YValidRML200On=allRML200On(~isnan(allRML200On));
allRMValidOn=nan(4000,2);
allRMValidOn(1:size(XValidRMControlOn),2)=XValidRMControlOn;
allRMValidOn(1:size(YValidRML200On),1)=YValidRML200On;
HistRMAllOn=hist3(allRMValidOn,[30 30]);
figure(126);
contour(HistRMAllOn,'fill','on');
colorbar; colormap(jet(40));
xlabel('Control');ylabel('Laser');
title('On Reponse');
figure(127); surf(hist3(allRMValidOn,[30 30]),'FaceAlpha',1);view([115,54])
colorbar; colormap(jet(20));
xlabel('Control');ylabel('Laser');
title('Reducing Median Clusters On Reponse');
colorbar; colormap(jet(20));
saveas(figure(127),[pwd '/New figures Jesus/ControlandL200/Surf ON FilterByMedian_ControlsVSlaser200.emf']);
saveas(figure(127),[pwd '/New figures Jesus/ControlandL200/Surf ON FilterByMedian_ControlsVSlaser200.pdf']);
saveas(figure(126),[pwd '/New figures Jesus/ControlandL200/Contour ON FilterByMedian_ControlsVSlaser200.emf']);
saveas(figure(126),[pwd '/New figures Jesus/ControlandL200/Contour ON FilterByMedian_ControlsVSlaser200.pdf']);

%% off
ReducingClustersRMControlOff=RspControlOff(:,SmallerMedianClusters);
allRMControlOff=ReducingClustersRMControlOff(:);
ReducingClustersRML200Off=RspL200Off(:,SmallerMedianClusters);
allRML200Off=ReducingClustersRML200Off(:);
boxplot([allRMControlOff allRML200Off],OffLabels)
title(["Off Response p=" ,signrank(allRMControlOff,allRML200Off)])
XValidRMControlOff=allRMControlOff(~isnan(allRMControlOff));
YValidRML200Off=allRML200Off(~isnan(allRML200Off));
allRMValidOff=nan(4000,2);
allRMValidOff(1:size(XValidRMControlOff),2)=XValidRMControlOff;
allRMValidOff(1:size(YValidRML200Off),1)=YValidRML200Off;
HistRMAllOff=hist3(allRMValidOff,[30 30]);
figure(125);
contour(HistRMAllOff,'fill','on');
colorbar; colormap(jet(20));
xlabel('Control');ylabel('Laser');
title('Off Reponse');
figure(124); surf(hist3(allRMValidOff,[60 60]),'FaceAlpha',1);view([115,54])
colorbar; colormap(jet(20));
xlabel('Control');ylabel('Laser');
title('Reducing Median Clusters Off Reponse');
colorbar; colormap(jet(20));
saveas(figure(124),[pwd '/New figures Jesus/ControlandL200/Surf OFF FilterByMedian_ControlsVSlaser200.emf']);
saveas(figure(124),[pwd '/New figures Jesus/ControlandL200/Surf OFF FilterByMedian_ControlsVSlaser200.pdf']);
saveas(figure(125),[pwd '/New figures Jesus/ControlandL200/Contour OFF FilterByMedian_ControlsVSlaser200.emf']);
saveas(figure(125),[pwd '/New figures Jesus/ControlandL200/Contour OFF FilterByMedian_ControlsVSlaser200.pdf']);
%% pie chart proportion
figure(115),pie3([sum(SmallerMedianClusters==1),sum(SmallerMedianClusters==0)]);
legend('reducing','increasing')
saveas(figure(115),[pwd '/New figures Jesus/all/Pie Chart reducing increasing.emf']);
saveas(figure(115),[pwd '/New figures Jesus/all/Pie Chart reducing Increasing.pdf']);

YControlRMstdOn=nanstd(RspControlOn(:,SmallerMedianClusters));
yL200RMstdOn=nanstd(RspL200On(:,SmallerMedianClusters));
%%boxplot([yControlRMstdOn yL200RMstdOn]);






