% Script to convert .mat files containing a single tetrode recording from
% Jurij and Adriano to .bin.
clearvars
%% File selection
dataFolder = 'C:\Users\caste\OneDrive\Desktop\Jurij Matlabscripts\tetrode8';
dataFileName = 'MEL12_PG9_REM_LOW_RESP_Tetrode_8.mat';
load(fullfile(dataFolder,dataFileName),'EEG')
%% Data scaling
% The data need to be transformed from any current scale into the 16 bit
% depth integers originated by the ADC. Normally, the ADCs use a 16 bit
% representation (2^16). After inspecting in EEGProcessing, you might get a
% slight idea of what the current scaling is. For what Melina and I,
% Emilio, saw, the scaling is from -50 mV to 50 mV with the data in uV,
% which gives us 100 mV -> 100.000 uV (1*10^5)
m = (2^16)/1e4; 

channelSelection = [1,3,5,7];
data = EEG.Data;
data = data(channelSelection,:);
dataScaled = int16(data.*m);

%% File writting
% To write the data into a binary file, we need to first define the
% location and name of the output file.
outputDirectory = 'F:\Kilosorting\MB-19031_PAC Genia\PACtetrode8';
outputFileName = 'MEL12_PG9_REM_LOW_RESP_Tetrode_8.bin';

fid = fopen(fullfile(outputDirectory,outputFileName),'w');
fwrite(fid,dataScaled,'int16');
fclose(fid);

%% Running kilosort
% First of all, you need to define some initial variables. For example, the
% location of your data which is in the 'outputDirectory' variable in the
% previous section. You need to have the configuration file together with
% the channel map!!
pathToYourConfigFile = outputDirectory; % take from Github folder and put it somewhere else (together with the master_file)
run(fullfile(pathToYourConfigFile, 'configFileSingleTetrode.m'))
rootH = pathToYourConfigFile;
ops.fproc       = fullfile(rootH, 'temp_wh.dat'); % proc file on a fast SSD
% ops.chanMap = fullfile(pathToYourConfigFile, 'neuropixPhase3A_kilosortChanMap.mat');

ops.trange = [0 Inf]; % time range to sort
ops.NchanTOT    = 4; % total number of channels in your recording

% the binary file is in this folder
rootZ = pathToYourConfigFile;


fprintf('Looking for data inside %s \n', rootZ)

% is there a channel map file in this folder?
fs = dir(fullfile(rootZ, 'chan*.mat'));
if ~isempty(fs)
    ops.chanMap = fullfile(rootZ, fs(1).name);
end

% find the binary file
fs          = [dir(fullfile(rootZ, '*.bin')) dir(fullfile(rootZ, '*.dat'))];
ops.fbinary = fullfile(rootZ, fs(1).name);

% preprocess data to create temp_wh.dat
rez = preprocessDataSub(ops);

% time-reordering as a function of drift
rez = clusterSingleBatches(rez);
save(fullfile(rootZ, 'rez.mat'), 'rez', '-v7.3');

% main tracking and template matching algorithm
rez = learnAndSolve8b(rez);

% final merges
rez = find_merges(rez, 1);

% final splits by SVD
rez = splitAllClusters(rez, 1);

% final splits by amplitudes
rez = splitAllClusters(rez, 0);

% decide on cutoff
rez = set_cutoff(rez);

fprintf('found %d good units \n', sum(rez.good>0))

% write to Phy
fprintf('Saving results to Phy  \n')
rezToPhy(rez, rootZ);

% discard features in final rez file (too slow to save)
rez.cProj = [];
rez.cProjPC = [];

% save final results as rez2
fprintf('Saving final results in rez2  \n')
fname = fullfile(rootZ, 'rez2.mat');
save(fname, 'rez', '-v7.3');

