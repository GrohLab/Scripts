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
    return
end
% If RollerSpeed.mat doesn't exist in the experiment folder, then create
% it.
rsfName = fullfile(dataDir, "RollerSpeed" + string(expDate) + ".mat");
if ~exist(rsfName, 'file')
    save(rsfName, "vf", "rollTx", "rollFs", "atTimes", "itTimes", "rx")
end
%% Figure time! 
fFig = figure(); 
plot(rollTx(1:end-1) - rollTx(1), vf, "DisplayName", "Roller velocity (speed)")
hold on; stem(itTimes(:,1), ones(size(itTimes, 1), 1), "DisplayName", "Intan triggers")
stem(atTimes, -ones(size(atTimes, 1), 1), "DisplayName", "Arduino triggers")
lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best")
title("Roller velocity and triggers"); xlabel("Time [s]"); 
ylabel("Roller velocity [encoder step per second]")
set(gca, "Box", "off", "Color", "none")
saveFigure(fFig, fullfile(dataDir, "Roller velocity and triggers"), 1);
%% Trigger average of the roller speed.
timeLapse = [-1, 2];
[~, vStack] = getStacks(false, atTimes * rollFs, 'on', timeLapse, rollFs,...
    rollFs, [], vf); vStack = squeeze(vStack);
plotEEGchannels(vStack', [], diff(timeLapse), rollFs, 1, -timeLapse(1))

[Nts, NTa] = size(vStack);
[ms, bs] = lineariz([1, Nts], timeLapse(2), timeLapse(1));
% Stack time axis --> st T x
stTx = (1:Nts)'*ms + bs;
% Spontaneous logical flag.
spontFlag = stTx < 0;
% No previous movement before nor after the stimulus
rmsTh = 0.85;
% Logical Flag to exclude trials with 'too much movement' before or after
% the trigger.
excludeFlag = rms(vStack(spontFlag,:)) > rmsTh | rms(vStack) > 0.9;
% Previous movement before 
% excludeFlag = rms(vStack(spontFlag,:)) < 0.85;

delayFlags = true(NTa,1); Nccond = 1;
% Number of triggers (Number of alignment points) Na excluding the
% indicated trials.
Na = sum(delayFlags & ~excludeFlag(:));

% Loop to get the maximum value of the roller velocity.
rngRollSpeed = cell(Nccond,1);
responseWindow = [0, 0.4];
responseFlags = stTx >= responseWindow(1) & stTx <= responseWindow(2);
for ccond = 1:Nccond
    rngRollSpeed{ccond} =...
        max(abs(vStack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end
% Maximum number of trials
maxNt = max(cellfun(@numel, rngRollSpeed));
% Placing the trial's maximum speed in a matrix
speedsMat = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    rngRollSpeed, fnOpts{:}); speedsMat = cat(2, speedsMat{:});
speedsMat = speedsMat * en2cm;
% Cutting on different thresholds for all trials
thetaSpeed = 0.1:0.1:3;
moveFlag = speedsMat > thetaSpeed;
probMove = sum(moveFlag)./Na; probMove = squeeze(probMove);
probFig = figure; plot(thetaSpeed, probMove'); ylim([0,1])
set(gca, "Box", "off", "Color", "none")
generalProb = nnz(moveFlag)./numel(moveFlag);
puffIntensity = strsplit(string(dataDir), "\");
puffIntensity = puffIntensity(end);
lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
ylabel("Trial proportion / Movement probability")
xlabel("Speed threshold \theta [cm/s]")
title(sprintf("Trial proportion with elicited movement\n(General prob: %.2f)",...
    generalProb))
%% Save?
sveFigAns = questdlg("Save figure?", "Save", "Yes", "No", "Yes");
probFigName = ...
    sprintf('Movement probability EX%.2f RW%.1f - %.1f s',...
    rmsTh, responseWindow);
if strcmpi(sveFigAns, "Yes")
    saveFigure(probFig, fullfile(dataDir, probFigName), 1)
end