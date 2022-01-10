% Reading the roller positions in a matrix of Nx2 (rp) where the first
% column is the position and the second the time in microseconds, a cell
% array with the trigger times (tTimes), and the corrected values from the
% arduino rollTrigTable
% Microseconds
us = 1e-6;
% cell- and arrayfun auxiliary variable.
fnOpts = {'UniformOutput', false};
% Encoder to centimeters per second.
en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*rollFs;
%% Choosing the file
rfName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211207\1.6bar\Roller_position2021-12-07T18_42_05.csv";
% rfName = uigetdir("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\",...
%     "Choose directory to work with");
[rp, tTimes, rollTrigTable] = readRollerPositionsFile(rfName);
%% Reading and preprocessing roller speed
% Get the folder and file name of the experiment to get the video frame
% rate to re-sample the position signal.
[dataDir, rollBaseName] = fileparts(rfName);
expDate = extractAfter(rollBaseName, "Roller_position");
% Split in date and time get the video 10 seconds ahead
dateFormStr = 'yyyy-MM-dd''T''HH_mm_ss';
expDate = datetime(expDate, 'InputFormat', dateFormStr);
vidDate = expDate + seconds(10); vidDate.Format = dateFormStr;
vfName = fullfile(dataDir, "roller" + string(vidDate) + ".avi");
vidObj = VideoReader(vfName); fr = vidObj.FrameRate;
% Computing a regular-interval time array in seconds to interpolate in
% between the asynchronous roller position values.
rollFs = fr; rollTx = (rp(1,2)*us:1/rollFs:rp(end,2)*us)';
rx = interp1(rp(:,2)*us, rp(:,1), rollTx, 'pchip');
%TODO: Beautify this figure
figure; plot(rp(:,2)*us, rp(:,1), rollTx, rx);
title("Roller position")
% Now, we compute the roller speed. Finally. We get the speed bc its easier
% to locate the small movements of the mouse.
v = diff(rx); 
[b, a] = butter(10, (2*18)/rollFs, "low");
vf = filtfilt(b, a, v);
%TODO: Beautify figure with butter and flies, or fries.
figure; plot(rollTx(1:end-1), [v, vf]);
%% Getting trigger times
% The point of this section is to equalize the triggers from the arduino
% with the trigger signals recorded from the Intan board.
expDate.Format = dateFormStr; atTimes = tTimes{:}.*us;
% Similar strategy as with the video but without addin 10 seconds.
fID = fopen(fullfile(dataDir,...
    "TriggerSignals" + string(expDate) + ".bin"), "r");
trig = fread(fID, [2, Inf], 'uint16=>int32'); [~] = fclose(fID);
trig = median(trig, 2) - trig; %trig = int16(trig);
figure; plot(trig'); 
% Sampling frequency of the intan board: 30 kHz. We also extract the
% trigger onsets of each pulse using the StepWaveform object.
fs = 3e4; pObj = StepWaveform(trig(1,:), fs);
pSubs = pObj.subTriggers; itTimes = pSubs./fs;
% Moving all arduino times from their arbitrary time back to zero.
atTimes = atTimes - rollTx(1);
% Removing repeated values from the trigger times.
atTimes([false;diff(atTimes(:,1)) < 1],:) = [];
% If the arduino trigger times are different from the intan, we need to
% fix it.
if size(itTimes,1) ~= size(atTimes)
    % Oh no... they are different. we can exclude the intan extra times.
    % But if you really feel up for it, you can estimate the time of
    % arduino. I can search for it in my MATLAB history.
    fprintf(1, 'Hm... What should we do?');
end
% If RollerSpeed.mat doesn't exist in the experiment folder, then create
% it.
rsfName = fullfile(dataDir, "RollerSpeed" + string(expDate) + ".mat");
if ~exist(rsfName, 'file')
    save(rsfName, "vf", "rollTx", "rollFs", "atTimes", "itTimes", "rx")
end
%% Trigger average of the roller speed.
timeSpan = [-1, 2];
[~, vStack] = getStacks(false, atTimes * rollFs, 'on', timeSpan, rollFs,...
    rollFs, [], {vf});