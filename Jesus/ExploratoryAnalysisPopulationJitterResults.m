
%% analyze population jitter results practice
clear all; close all; % clear everything
load ('D:\190702_Jittering_3720_1520_1520\jesus_emilio_jittering1_PopulationSpikeAnalysis.mat'); % load files


%%
% PopRelativeSpikeTimes{x}{y,z}
% x= condition number, y=cluster number, z = trial numbers
fig=figure
n=7    %
H=[]
for i=1:size(PopRelativeSpikeTimes{n},1)  %for all neurons
    N=size(PopRelativeSpikeTimes{n}(i,:),2) %number of trials
    spikes=cell2mat(PopRelativeSpikeTimes{n}(i,:));
    binsize=1;%in ms
    hbins=(-timeLapse(1)*1000):binsize:(timeLapse(2)*1000); %in ms
    h=histogram(spikes*1000,hbins);
    figure(fig)
    H(i,:)=h.Values;
    
end
%%

figure
plot(hbins(2:end)',H)  
xlim([-100 150])


%% step through all of the units
figure
for i=1:63
    plot(hbins(2:end)',H(i,:))
    title(i)
    pause
end
xlim([-100 150])
%% how many response to whisker?

%% how many reponses to light? 



%% split ex and inhib by wavefrom, independent of spike time data
%first would be ideal to split spikes depending on first component
% waveform duration and ISI. as it has been related putative type 1 neurons
% wider first component. and putative type 2 with shorter first component
% duration and also with higher spike frequency.



%% load the data saved in the other script

%% basic plotting of psth by condition
%driven by light: Psth triggered by Laser stimulation.

%driven by whisker: PSTH triggered by piezo in...
%at what time lag: Control, 1ms, 10ms, 50ms,100ms,200ms.

%sequence of analysis steps to be done for each condition

%% whisker affected by light
% how to do this?
% spike count in a defined window
% response probability
% characterization of spike distribution
% what do we want to measure?
% KS
% what tests are appropriate for comparing distribution

%% spontaneous vs. stimulus driven

% for spontaneous activity
% i am not really sure what information spontaneous activity will give us in VPM.
%  mean frequency, or instanteous frequency (taking care of After hiperpolarization currents
% that can be observerd after Piezo.

%for stimulus driven responses

%response probability ... in each condition.
% number of spikes per cluster in each bin of the PSTH ( being the bin of  1ms)
% response probability of population responsive clusters
% Population responsive clusters PSTH: number of spikes per bin (bin 0.5ms?)

%% first spike statistics

%standard deviation of first spike, cumulative fraction of first spike,
% using the first spike: reverse spike average
% probability of spikes in the first bin/total number of spikes

%% all spike statistics...
% distribution comparation of cumulative fraction.
% number of total spikes in each condition
%
%cumulative fraction all spikes in 50ms windows after piezo onset.
%width of Psth in each condition (min - max latency, depicted in box and
%whiskers plot)
% psth could be splitted in two components: first component 0-20/25ms, second
% component 25-50/60ms (it is has been related with AMPA and NMDA
% components) and compare in different laser condictions and control.


%% bursting
% interval spike histogram for each cluster. temporal plot to check burst
% appearance. Number of spikes in each burst.



