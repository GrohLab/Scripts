% Reading the roller positions in a matrix of Nx2 (rp) where the first
% column is the position and the second the time in microseconds, a cell
% array with the trigger times (tTimes), and the corrected values from the
% arduino rollTrigTable
us = 1e-6;
rfName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211209\Roller_position2021-12-09T18_35_15.csv";
[rp, tTimes, rollTrigTable] = readRollerPositionsFile(rfName);
%% Roller speed
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
vf = filtfilt(b, a, rs);
%TODO: Beautify figure with butter and flies.
figure; plot(rollTx(1:end-1), [v, vf]);
%% Trigger times
% The point of this section is to equalize the triggers from the arduino
% with the trigger signals recorded from the Intan board.
