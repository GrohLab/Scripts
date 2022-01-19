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
rfName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211206\1.8bar\Roller_position2021-12-06T18_46_05.csv";
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
%% Trigger average of the roller speed.
timeLapse = [-1, 2];
responseWindow = [0, 0.4]; % will stay the same right? Yes

[~, vStack] = getStacks(false, atTimes * rollFs, 'on', timeLapse, rollFs,...
    rollFs, [], vf);

[~, Nts, NTa] = size(vStack);
[ms, bs] = lineariz([1, Nts], timeLapse(2), timeLapse(1));

% Stack time axis --> st T x
stTx = (1:Nts)'*ms + bs;

rmsTh1 = 0.85;
rmsTh2 = 0.9;
thetaSpeed = 0.1:0.1:3;
[genProbSpeed, probMatSpeed] = getGeneralProb(vStack, thetaSpeed, en2cm, ...
    rmsTh1, rmsTh2 );

%Plot EEG channels
plotEEGchannels(vStack', [], diff(timeLapse), rollFs, 1, -timeLapse(1))

%Plot General Probability
probFig = figure; 
plot(thetaSpeed, probMatSpeed'); 
ylim([0,1])
set(gca, "Box", "off", "Color", "none")
%generalProb = nnz(moveFlag)./numel(moveFlag);
puffIntensity = strsplit(string(dataDir), "\");
puffIntensity = puffIntensity(end);
lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
ylabel("Trial proportion / Movement probability")
xlabel("Speed threshold \theta [cm/s]")
title(sprintf("Trial proportion with elicited movement\n(General prob: %.2f)",...
    genProbSpeed))
%% Save?
sveFigAns = questdlg("Save figure?", "Save", "Yes", "No", "Yes");
probFigName = ...
    sprintf('Movement probability EX%.2f RW%.1f - %.1f s',...
    rmsTh1, responseWindow);
if strcmpi(sveFigAns, "Yes")
    saveFigure(probFig, fullfile(dataDir, probFigName), 1)
end

% If RollerSpeed.mat doesn't exist in the experiment folder, then create
% it.
rsfName = fullfile(dataDir, "RollerSpeed" + string(expDate) + ".mat");
if ~exist(rsfName, 'file')
    save(rsfName, "vf", "rollTx", "rollFs", "atTimes", "itTimes", "rx", "genProbSpeed")
end
%% DeepLabCuts - whisker movements as a reaction. 
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
sNames = ["LeftWhiskers","RightWhiskers","Nose"];
[~, dlcStack] = getStacks(false, round(itTimes * fr), 'on', timeLapse,...
    fr, fr, [], lw, rw, nose);

lwstack = squeeze(dlcStack(1,:,:));
rwstack = squeeze(dlcStack(2,:,:));
nosestack = squeeze(dlcStack(3,:,:));

thetaAnglWh = 0.5:0.5:40;
thetaAnglNose = 0.2:0.2:10;
rmsTh3 = 0.85;
rmsTh4 = 6;

[genProbLW, probMatLW] = getGeneralProb(lwstack, thetaAnglWh, 1, ...
    rmsTh3, rmsTh4);
[genProbRW, probMatRW] = getGeneralProb(rwstack, thetaAnglWh, 1, ...
    rmsTh3, rmsTh4);
[genProbNose, probMatNose] = getGeneralProb(nosestack, thetaAnglNose, 1, ...
    rmsTh3, rmsTh4);

%% Loop to get the maximum angle of the different bodyparts.

angllw = cell(Nccond,1);
anglrw = cell(Nccond,1);
anglnose = cell(Nccond,1);
% responseWindow = [0, 0.4];
% just takes the samples from 0 to 0.4
% responseFlags = stTx >= responseWindow(1) & stTx <= responseWindow(2);

for ccond = 1:Nccond
  angllw{ccond} =...
         max(abs(lwstack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end

for ccond = 1:Nccond
  anglrw{ccond} =...
         max(abs(rwstack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end

for ccond = 1:Nccond
  anglnose{ccond} =...
         max(abs(nosestack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end

% Maximum number of trials
maxNt = max(cellfun(@numel, angllw));
% Placing the trial's maximum angle in a matrix
anglMatLW = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    angllw, fnOpts{:});
anglMatRW = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    anglrw, fnOpts{:});
anglMatNose = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    anglnose, fnOpts{:});

anglMatLW = cat(2, anglMatLW{:});
anglMatRW = cat(2, anglMatRW{:});
anglMatNose = cat(2, anglMatNose{:});

anglMatLW = anglMatLW * en2cm;
anglMatRW = anglMatRW * en2cm;
anglMatNose = anglMatNose * en2cm;

% Cutting on different thresholds for all trials

anglFlagLW = anglMatLW > thetaAnglWh;
anglFlagRW = anglMatRW > thetaAnglWh;
anglFlagNose = anglMatNose > thetaAnglNose;

probAnglLW = sum(anglFlagLW)./Na;
probAnglRW = sum(anglFlagRW)./Na;
probAnglNose = sum(anglFlagNose)./Na;

probAnglLW = squeeze(probAnglLW);
probAnglRW = squeeze(probAnglRW);
probAnglNose = squeeze(probAnglNose);

generalAnglProbLW = nnz(anglFlagLW)./numel(anglFlagLW);
generalAnglProbRW = nnz(anglFlagRW)./numel(anglFlagRW);
generalAnglProbNose = nnz(anglFlagNose)./numel(anglFlagNose);

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
    
    figname = ['SumOfSignals', num2str(sNames(cbp)), '.emf'];
    saveFigure(efig(cbp), fullfile(dataDir, figname), 1);
end