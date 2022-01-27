% HBIGs 2022 in vivo electrophysiology analysis tutorial

% Welcome to the HBIGS tutorial. In this tutorial will be be covering the
% basics of the kinds of analysis we do in the Groh Lab, specifically
% multi-unit extracellular electrophysiology.

% Whilst the script should run in full, we recommend (especially for those
% unfamiliar with MATLAB or programming) to run the script section by
% section to really get a feel for how and why each step is necessary for answering the questions that we're interested in. 
% You can do this run each section by either clicking th 'Run Section' icon in the
% Editor section of the taskbar, or by pressing Strg + Enter.
% To run a section and then automatically advance to the next section in
% the script there is also a 'Run and Advance' button, or alternatively
% Strg + Shift + Enter will do this.



%% Firstly, we need a data directory to work in:

dataDir = 'Z:\HBIGS';

% We can see what files are in the current directory with the ls functino:

ls(dataDir)

%% Concatenating and converting smrx files to bin file.

% The .smrx files contain the data that was collected during the experiment.
% Therefore we need to access both the voltage traces from the probes, and
% therefore the traces of any non-neural paramters we may be interested in
% (e.g. TTL traces, respiration rate, etc.)


% To take votlage traces from multiple electrodes and and assign the
% spikes from these into respective neurons, we need to put these traces
% through a 'Spike-sorter'; an algorithm that clusters similar spikes as
% belonging to these putative neurons (referred to as units). But first,
% we need to concatenate the separate smrx files together and convertthem into a single binary file to be read by the
% spike-sorting software.

 
 msmrx2bin('Z:\HBIGS','HBIGS_Tutorial')
 
 % The chronological order of the recording are denoted by 'P', where 'P1_'
 % denotes the first recording of the experiment.
 %% Now that we have the bin file, we can run our spike-sorter. 
 
 
 cd C:\Users\NeuroNetz\Documents\GitHub\Kilosort2_5
 
 kilosort
 
 % probe file = 'Z:\Ross\ProbeFiles\Corrected_H3_ChanMap.mat'
 % sampling frequency = 3.003003003003003e+04
 