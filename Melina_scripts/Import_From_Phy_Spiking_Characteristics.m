% import data from Phy and perform very basic spiking characterizing and
% save the data for later plotting %r.mease and m.castelanelli
clear all
%give a directory where the Phy curated data lives
startpath='F:\Kilosorting\'
datadir = uigetdir(startpath,  'Select a data directory')

%or hardcode
%datadir = 'F:\Kilosorting\#7-2019_G29-16\Recording-1';



%%
outputName='phySpikeOutput';
removeNoise=true %get good, mua, noise, but only want to use good sorted spikes for further analysis
%import the data
importPhyFiles(datadir,outputName,removeNoise)
cd(datadir)
%%
%spikes=sortedData{1,2}  % pick out the spike times from the first ID (first row, second column) but
%want all good spike times
matlabfile = [outputName '_all_channels'];
load(matlabfile);

%from when importPhyFiles did not kick out bad clusters
%badsIdx = cellfun(@(x) x==3,sortedData(:,3)); %old command to remove bads
%Spikes=sortedData(find(~badsIdx),2); %Spikes means spike times fromm all "good" (=[1]) spikes
%ClusterIDs=sortedData(find(~badsIdx),1);% ID´s of all "good" spikes

Spikes=sortedData(:,2);
ClusterIDs=sortedData(:,1)

%%
%Calculating the spike rates (=how often the neuron fires)
Rates=zeros(size(Spikes));
Isis={};
%do rate and ISI calculations per cluster
for i=1:numel(Spikes)
    spikes=Spikes{i};
    T=max(spikes)  %total time of experiment (needs to be found out somehow) This is cheating currently, just took the last spike as recording length!
    Rates(i)= length(spikes)/T   %calculate rate. length means the length of a vector (amount of events. e.g. 13 -> means neurons fire 13 times per second thus 13 Hz
    Isis{i}=diff(spikes)*1000;%looking at ms. ISI is difference between spike times
end

save SpikingParameters Spikes ClusterIDs Rates Isis datadir
