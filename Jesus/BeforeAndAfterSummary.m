% do connected categorical scatter plots

% add titles and labels
% add variables to change the desired window
% measure median for each unit
%add the box plot
%add experimental id
%think about statistical test to see effect on single neuron level.
%next time with rebecca add dkl
%need to find the appropriate nonparametric 2 factor anova



rsp=ClustersControlOnset*1000; %now in ms
%find all unwanted spikes
exclude=[find(rsp<2);find(rsp>20)];
rsp(exclude)=nan;
y=nanmedian(rsp);


rsp=ClustersL200Onset*1000; %now in ms
%find all unwanted spikes
exclude=[find(rsp<2);find(rsp>20)];
rsp(exclude)=nan;
y1=nanmedian(rsp);

x=ones(size(y))*1.2;
x1=ones(size(y1))*1.8;

figure
hold on
for i=1:numel(x1)
   plot([1 2],[y(i) y1(i)],'-o') 
end
xlim([0 3])

title(('test')%signrank(y,y1))  %for paired data
 

%%
% measure median for each unit
rsp=ClustersControlOnset*1000; %now in ms
%find all unwanted spikes
exclude=[find(rsp<2);find(rsp>20)];
rsp(exclude)=nan;
y=nanstd(rsp);

rsp=ClustersL200Onset*1000; %now in ms
%find all unwanted spikes
exclude=[find(rsp<2);find(rsp>20)];
rsp(exclude)=nan;
y1=nanstd(rsp);

x=ones(size(y));
x1=ones(size(y1))*2;

figure
hold on
for i=1:numel(x)
   plot([1 2],[y(i) y1(i)],'-o') 
end
xlim([0 3])

title(signrank(y,y1))  %for paired data
 

%

%combine with boxplot

%measure DKL for each unit (before after)
%violin plots