% Melina script for data loading and notch filtering
clear all  %start with empty workspace

%Loading the recording; this also defines the filename using the file
%manually selected
clear all
rawDir='S:\Rebecca\RawData'
saveDirRoot='S:\Rebecca\AnalyzedData';

cd(rawDir)
read_Intan_RHD2000_file  %to do is to remove gui and just do a loop to import all data in directory
relativePath=path((numel(rawDir)+2):end);
savePath=[saveDirRoot '\' relativePath];
cd(saveDirRoot)
if ~isdir(relativePath)
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(savePath,relativePath);
    if SUCCESS,
        display(['made new directory ' savePath])
    else display(MESSAGE)
    end
else
    display(['directory ' savePath ' already exists'])
end
cd(savePath)

%% then filter out 50 Hz line noise and filter for spikes
fs = frequency_parameters.amplifier_sample_rate; %original rate
[rows, cols] = size(amplifier_data);
if cols < rows
    amplifier_data = amplifier_data';
    aux = rows;
    rows = cols;
    cols = aux;
end

dataCell = cell(rows,1);dataCell_sp = cell(rows,1); %empty structs for filtered data
ds_factor=10; %by what factor should signal be downsampled by? 
fs_ds=fs/ds_factor;
cutFreq = [350, 1400]; %is this too narrow?
for cch = 1:min(size(amplifier_data))
    fprintf(1,'Dealing with channel %d\n', cch)
    ds_signal=resample(amplifier_data(cch,:),1,ds_factor);
    dataCell(cch) = {iir50NotchFilter(ds_signal,fs_ds)};
    dataCell_sp(cch) = {iirSpikeFilter(dataCell{cch},fs_ds,cutFreq)};
   % dataCell_sp(cch) = {iirSpikeFilter(ds_signal,fs_ds,cutFreq)};
   
end
cd(savePath)


%%
% we need to have a general strategy for reading in the commands NEEDS work
laserSignal = board_dig_in_data(4,:);  %% keep the laser command...
clearvars -except laser* fs cols rows dataCell dataCell_sp matfilename path savePath
save([matfilename '_denoise'],'dataCell'); %raw data
save([matfilename '_highpass'],'dataCell_sp'); %high pass for spikes
analysisFile = [matfilename 'analysis']; %just the path and name as string
%THIS WILL OVERWRITE ANY ANALYSIS THAT HAS BEEN DONE, NEEDS TO BE FIXED
save(analysisFile,'matfilename', 'path','fs','analysisFile','laserSignal','savePath') %initial save of analysis file to store
