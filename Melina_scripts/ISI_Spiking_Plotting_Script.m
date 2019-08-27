%% load in spiking data and precalculated rates etc...
clear all
[file, path, filterindex] = ...
    uigetfile('*.mat', '´Select a file with spike rates and ISIs', 'MultiSelect', 'off');
load( [path,file])


%% interspike intervals for one ID (one neuron)
close all
%Make this into a loop for all clusters
isis=Isis{1}
close all
%make a loop to make isis histograms
histogram(log10(isis),200)
xtick_ms=10.^get(gca,'XTick');
set(gca,'XTickLabel',xtick_ms)
title(ClusterIDs{1})
xlabel('ISI in ms')
ylabel('firing events')

%% interspike intervals for all clusters (using a loop)
close all

for i=1:numel(Isis)
    figure
    isis=Isis{i};
    histogram(log10(isis),200)
    xtick_ms=10.^get(gca,'XTick');
    set(gca,'XTickLabel',xtick_ms)
    title([ 'ClusterId=' ClusterIDs{i}  ' Cluster number=' num2str(i)])
    xlabel('ISI in ms')
    ylabel('events')
    
end

%% exercise-- make a histogram out of spike rates (Rates)
figure

histogram(Rates,20)  % THIS LOOKS SUPER STRANGE!! he can not plot all in one figure

%  title([ 'ClusterId=' ClusterIDs{i}  ' Cluster number=' num2str(i)])
xlabel('spike rate in Hz')
ylabel('count') %count means number of contributing neurons in that Hz

ylim([0 10])


%% just look at the spike times for one/all cluster
%Visualize spikes over time --> raster plot for a given chunk of time

close all
minTime=100 %in seconds
maxTime=200  %in seconds  %checkt that this is within the range of the recording

%spikes=Spikes{1}; %grab spike times from first cluster
figure
level=0
for i=1:numel(Spikes)  %grab spike times from all "good" cluster
    level=level+1;
    spikes=Spikes{i};
    spikes=spikes(spikes>minTime & spikes<maxTime);
    y=ones(size(spikes))*level;
    plot(spikes,y,'LineStyle','none','Marker','.','Markersize',10)
    title('spike times per cluster')
    xlabel('time in sec')
    ylabel('clusters')
    hold on
end
ylim([0 i+1])   %add some space before first cluster and last cluster that it looks nicer in the histogram


%load whisker information if available
try 
    load WhiskerDataRough
    indices=find(t_dig>minTime & t_dig<maxTime);
    plot(t_dig(indices), puff(indices)*( i+1),'k')
catch
end



%% set up look to look at cross correlation between spike trains
    
%% time binned rate for one cluster

spikes=Spikes{1};   %we use here spike times from cluster 1
figure
binsize=1 %in sec (timewindow which groups spike firing of neuron 1
%H=histogram(spikes,'BinWidth',binsize)%,making histogram with variable width (histogram of spike times cluster 1 grouped in bins) 'BinWidth',2 !!!
H=histogram(spikes,'BinWidth',binsize, 'Normalization','countdensity')
title('time binned rates #1')
xlabel('time in sec')
ylabel('spike rate in Hz') %spike counts normalized and bin size is in seconds
R=H.Values;

%% time binned rate for each cluster

for    i=1:numel(Spikes)
    spikes=Spikes{i};
    figure
    binsize=1 %in sec
    H=histogram(spikes,'BinWidth',binsize, 'Normalization','countdensity')
  
    title([ 'ClusterId=' ClusterIDs{i}  ' Cluster number=' num2str(i)])
    xlabel('time in sec')
    ylabel('spike rate in Hz') %spike counts normalized and bin size is in seconds
    R=H.Values;
    
end

%%EX figure out how to make line of plot thicker to make it easier to
        %see!!!


%%
close all
figure
%one compound plot of all rates
keepers=[];
binsize=5 %in sec

for i=1:numel(Spikes)   %for all the good neurons
% for i= [38 39 41 43 35 31 33 2 27 25 10 12]  % for specific selected
% good neurons  (tip in # not clusterID)
 
    if Rates(i)>.5
    spikes=Spikes{i};
    H=histogram(spikes,'BinWidth',binsize, 'Normalization','countdensity','DisplayStyle','stairs')
        %EX figure out how to make line of plot thicker to make it easier to
        %see!!!

    xlabel('time in sec')
    ylabel('spike rate in Hz') %spike counts normalized and bin size is in seconds
    R=H.Values;
    hold on
    keepers=[keepers ' ' ClusterIDs{i}];
    else
        display('Rate too low...')
    end
end
title(['Displaying clusters ' keepers])





  %% figure out how to keep the values of R for each neuron



%temporal window and then count all the spikes during this period of time
%%

%Overall spiking rate per unit -->  Number of spikes/time

