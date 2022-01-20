% Reading the roller positions in a matrix of Nx2 (rp) where the first
% column is the position and the second the time in microseconds, a cell
% array with the trigger times (tTimes), and the corrected values from the
% arduino rollTrigTable
% Microseconds
close all
us = 1e-6;
% cell- and arrayfun auxiliary variable.
fnOpts = {'UniformOutput', false};
%% Choosing the file
rfName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211206\0.2bar\Roller_position2021-12-06T17_53_22.csv";
% rfName = uigetdir("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\",...
%     "Choose directory to work with");

[rp, tTimes, rollTrigTable] = readRollerPositionsFile(rfName);
if length(tTimes) > 1
    tTimes = tTimes(2);
    fprintf(1, '"L" Trigger is omitted');
end
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
try
    vidObj = VideoReader(vfName);
catch
    vfName = dir(fullfile(dataDir, "*.avi"));
    if numel(vfName) > 1
        sel = listdlg("ListString", arrayfun(@(x) x.name, vfName,...
            'UniformOutput', false), "SelectionMode", "multiple");
        if isempty(sel)
            fprintf(1, 'Cancelling...\n')
            return
        end
        vfName = fullfile(vfName(sel).folder, vfName(sel).name);
    else
        vfName = fullfile(vfName.folder, vfName.name);
    end
end
fr = vidObj.FrameRate;
% Computing a regular-interval time array in seconds to interpolate in
% between the asynchronous roller position values.
rollFs = fr; 
rollTx = (rp(1,2)*us:1/rollFs:rp(end,2)*us)';
% Encoder to centimeters per second.
en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*rollFs;
% Roller time axis
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
% Sampling frequency of the intan board: 30 kHz. We also extract the
% trigger onsets of each pulse using the StepWaveform object.
fs = 3e4; 
pObj = StepWaveform(trig(1,:), fs);
pSubs = pObj.subTriggers; 
itTimes = pSubs./fs;
% Moving all arduino times from their arbitrary time back to zero.
atTimes = atTimes - rollTx(1);
% Removing repeated values from the trigger times.
atTimes([false;diff(atTimes(:,1)) < 1],:) = [];
% If the arduino trigger times are different from the intan, we need to
% fix it.
if size(itTimes,1) ~= size(atTimes,1)
    % Oh no... they are different. we can exclude the intan extra times.
    % But if you really feel up for it, you can estimate the time of
    % arduino. I can search for it in my MATLAB history.
    fprintf(1, 'Hm... What should we do?');
end
%% Figure time! 
figure; plot(trig'); 

fFig = figure(); 
plot(rollTx(1:end-1) - rollTx(1), vf, "DisplayName", "Roller velocity (speed)")
hold on; stem(itTimes(:,1), ones(size(itTimes, 1), 1), "DisplayName", "Intan triggers")
stem(atTimes, -ones(size(atTimes, 1), 1), "DisplayName", "Arduino triggers")
lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best")
title("Roller velocity and triggers"); xlabel("Time [s]"); 
ylabel("Roller velocity [encoder step per second]")
set(gca, "Box", "off", "Color", "none")
saveFigure(fFig, fullfile(dataDir, "Roller velocity and triggers"), 1);
%% Trigger average of behaviour signals.
timeLapse = [-1, 2];
responseWindow = [0, 0.4]; % will stay the same right? Yes

[~, vStack] = getStacks(false, atTimes * rollFs, 'on', timeLapse, rollFs,...
    rollFs, [], vf*en2cm);

%Plot EEG channels
plotEEGchannels(squeeze(vStack)', [], diff(timeLapse), rollFs, 1, -timeLapse(1))

[~, Nts, NTa] = size(vStack);
[ms, bs] = lineariz([1, Nts], timeLapse(2), timeLapse(1));
% Stack time axis --> st T x
stTx = (1:Nts)'*ms + bs;
spontFlag = stTx < 0;


% DeepLabCut part
[~, vfNameBase] = fileparts(vfName);
dlcffName = dir(fullfile(dataDir, string(vfNameBase)) + "*filtered.csv");
dlcffName = fullfile(dataDir, dlcffName.name);
dlcTable = readDLCData(dlcffName);
[a_bodyParts, refStruct] = getBehaviourSignals(dlcTable);
% Left whiskers
lw = mean(a_bodyParts{:,{'lw1', 'lw2', 'lw3', 'lw4'}},2);
lw = lw - mean(lw);
% Right whiskers
rw = mean(a_bodyParts{:,{'rw1', 'rw2', 'rw3', 'rw4'}},2);
rw = rw - mean(rw);
% Nose signal
nose = a_bodyParts{:,"nose"} - mean(a_bodyParts{:,"nose"});

% Trigger cut
sNames = ["LeftWhiskers", "RightWhiskers", "Nose", "RollerSpeed"];
[~, dlcStack] = getStacks(false, round(itTimes * fr), 'on', timeLapse,...
    fr, fr, [], lw, rw, nose);

puffIntensity = strsplit(string(dataDir), "\");
puffIntensity = puffIntensity(end);

%close all;
%% EEG-like plots for the nose and the two whisker-sets
efig = 1:size(dlcStack,1);

for cbp = 1:size(dlcStack,1)
    plotEEGchannels(squeeze(dlcStack(cbp,:,:))', [], diff(timeLapse),...
        fr, 1, -timeLapse(1));
    title(sNames(cbp))
    yticklabels(NTa:-1:1)
    efig(cbp) = figure();
    plot(stTx, squeeze(dlcStack(cbp,:,:)), "LineWidth", 0.5,...
        "Color", 0.75*ones(3,1)); hold on;
    plot(stTx, mean(dlcStack(cbp, :, :), 3), "LineWidth", 2, "Color", "k")
    title(sNames(cbp))
    lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
    
    figname = ['SumOfSignals', num2str(sNames(cbp))];
    saveFigure(efig(cbp), fullfile(dataDir, figname), 1);
    
end
%% Behaviour stack for calculating the probability once for all signals
behStack = dlcStack; behStack(4,:,:) = vStack;
% Exclude trials with spontaneous running, which will affect the face
% movements.
sSigma = squeeze(std(behStack(:, spontFlag, :), [], 2));
sigTh = [5; 5; 2; 2.5];
excludeFlag = sSigma > sigTh;
behStack = arrayfun(@(x) squeeze(behStack(x, :, :)), (1:size(behStack,1))',...
    fnOpts{:});
thSet = {0.5:0.5:40,... Left whiskers
    0.5:0.5:40,... Right whiskers
    0.2:0.2:10,... Nose
    0.1:0.1:3}; % Roller speed
Ns = length(behStack);
behResults = cell(Ns, 2);
genProb = zeros(Ns, 1);
for cs = 1:size(behStack,1)
    % Getting the maximum value for all signals within the response window
    mvpt = getMaxAbsPerTrial(behStack{cs}(:, ~excludeFlag(cs,:)),...
        responseWindow, stTx);
    % Compare the max value against a set of thresholds
    moveFlags = compareMaxWithThresh(mvpt, thSet(cs));
    behResults(cs,:) = {mvpt, moveFlags};
    % Getting the 'general probability' of movement (Area Under the Curve)
    genProb(cs) = getAUC(moveFlags);
end
%% Plotting the threshold progression
% probFigs = plotThetaProgress(moveFlags, thSet, sNames);
probFigs = arrayfun(@(x) plotThetaProgress(behResults(x, 2), thSet(x), ...
    sNames(x)), 1:Ns);
xlArray = [repmat("Angle [°]",1,3), "Speed [cm/s]"];
% xlArray = insertBefore(xlArray, "[", "\\theta ");
ttlString = sprintf("Trials crossing \\theta (%s)", puffIntensity);
arrayfun(@(x) xlabel(probFigs(x).CurrentAxes, xlArray(x)), 1:Ns)
arrayfun(@(x) title(probFigs(x).CurrentAxes, ttlString), 1:Ns)
arrayfun(@(x) saveFigure(probFigs(x), fullfile(dataDir, sNames(x) +...
    sprintf(" trial crossing (%s)", puffIntensity)), 1), 1:Ns)

%close all;
%% Save?
% If RollerSpeed.mat doesn't exist in the experiment folder, then create
% it.
rsfName = fullfile(dataDir, "BehaviourData" + string(expDate) + ".mat");
if ~exist(rsfName, 'file')
    save(rsfName, "vf", "rollTx", "rollFs", "atTimes", "itTimes", "rx",...
        "genProb", "thSet", "moveFlags", "mvpt", "behStack", "stTx")
end

%% EEG-like plots for the nose and the two whisker-sets
efig = 1:size(dlcStack,1);
for cbp = 1:size(dlcStack,1)
    plotEEGchannels(squeeze(dlcStack(cbp,:,:))', [], diff(timeLapse),...
        fr, 1, -timeLapse(1));
    title(sNames(cbp))
    efig(cbp) = figure(); 
    plot(stTx, squeeze(dlcStack(cbp,:,:)), "LineWidth", 0.5,...
        "Color", 0.75*ones(3,1)); hold on;
    plot(stTx, mean(dlcStack(cbp, :, :), 3), "LineWidth", 2, "Color", "k")
    title(sNames(cbp))
    lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
    
    figname = ['SumOfSignals', num2str(sNames(cbp)), ' filtered'];
    saveFigure(efig(cbp), fullfile(dataDir, figname), 1);
end

fprintf(1, "General probability: %.3f\n", genProb)
