function behaveAnalysis(dataDir)
% Reading the roller positions in a matrix of Nx2 (rp) where the first
% column is the position and the second the time in microseconds, a cell
% array with the trigger times (tTimes), and the corrected values from the
% arduino rollTrigTable
% Microseconds
us = 1e-6;
% cell- and arrayfun auxiliary variable.
fnOpts = {'UniformOutput', false};
figureDir = fullfile(dataDir, 'Figures');
axOpts = {"Box", "off", "Color", "none"};
mkdir(figureDir);
%% Roller speed and trigges from both arduino and intan
% Triggers
fullname = @(x) fullfile(x.folder, x.name);
trigFile = dir(fullfile(dataDir, 'Roller_position*.csv'));
expDate = getDates(trigFile, "Roller_position");
try
    [atTimes, ~, itTimes] = readAndCorrectArdTrigs(dataDir);
    atTimes = atTimes{1}; itTimes = itTimes{1};
catch
    trigFile = dir(fullfile(dataDir, "ArduinoTriggers*.mat"));
    load(fullname(trigFile));
    atTimes = atTimes{1}; itTimes = itTimes{1};
end
% Roller speed
try
    [~, vf, rollTx, fr] = createRollerSpeed(dataDir);
catch
    rsFiles = dir(fullfile(dataDir, "RollerSpeed*.mat"));
    load(fullname(rsFiles))
end
rollTx = rollTx{1};
en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fr;

%% Speed figure
fFig = figure(); 
plot(rollTx(1:end-1) - rollTx(1), vf, "DisplayName", "Roller velocity (speed)")
hold on; stem(itTimes(:,1), ones(size(itTimes, 1), 1), "DisplayName", "Intan triggers")
stem(atTimes, -ones(size(atTimes, 1), 1), "DisplayName", "Arduino triggers")
lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best")
title("Roller velocity and triggers"); xlabel("Time [s]"); 
ylabel("Roller velocity [encoder step per second]")
set(gca, axOpts{:})
saveFigure(fFig, fullfile(figureDir, "Roller velocity and triggers"), 1);
%% Parameters
timeLapse = [-1, 2];
responseWindow = [0, 0.4];
%% DeepLabCut signals
dlcffName = dir(fullfile(dataDir, "*filtered.csv"));
dlcffName = fullfile(dataDir, dlcffName.name);
dlcTable = readDLCData(dlcffName);
[a_bodyParts, ~] = getBehaviourSignals(dlcTable);

% Left whiskers
lw = mean(a_bodyParts{:,{'lw1', 'lw2', 'lw3', 'lw4'}},2);
lw = lw - mean(lw);
% Right whiskers
rw = mean(a_bodyParts{:,{'rw1', 'rw2', 'rw3', 'rw4'}},2);
rw = rw - mean(rw);
% Nose signal
nose = a_bodyParts{:,"nose"} - mean(a_bodyParts{:,"nose"});

% DLC Trigger cut
sNames = ["Ipsi-Whiskers", "Contra-Whiskers", "Nose", "RollerSpeed"];
[~, dlcStack] = getStacks(false, round(itTimes * fr), 'on', timeLapse,...
    fr, fr, [], lw, rw, nose);

% Roller speed trigger cut
[~, vStack] = getStacks(false, round(atTimes * fr), 'on', timeLapse, fr,...
    fr, [], vf*en2cm);

% Auxiliary variables
Nts = size(vStack,2);
[ms, bs] = lineariz([1, Nts], timeLapse(2), timeLapse(1));
% Stack time axis --> st T x
stTx = (1:Nts)'*ms + bs;
spontFlag = stTx < 0;


% Puff intensity from the folder name
puffIntensity = strsplit(string(dataDir), "\");
puffIntensity = puffIntensity(end);

%% Behaviour stack for calculating the probability once for all signals
behStack = dlcStack; behStack(4,:,:) = vStack;
% Exclude trials with spontaneous running, which will affect the face
% movements.
sSigma = squeeze(std(behStack(:, spontFlag, :), [], 2));
sigTh = [5; 5; 2; 2.5];
excludeFlag = sSigma > sigTh;
% This doesn't take into consideretion the nose
excludeFlag = any(excludeFlag(2.^[0:2],:),1); 
Net = sum(excludeFlag, 2);
behStack = arrayfun(@(x) squeeze(behStack(x, :, :)), (1:size(behStack,1))',...
    fnOpts{:});
thSet = {0.5:0.5:40,... Left whiskers
    0.5:0.5:40,... Right whiskers
    0.4:0.4:20,... Nose
    0.1:0.1:3}; % Roller speed
Ns = length(behStack);
behResults = cell(Ns, 2);
genProb = zeros(Ns, 1);
for cs = 1:size(behStack,1)
    % Getting the maximum value for all signals within the response window
    mvpt = getMaxAbsPerTrial(behStack{cs}(:, ~excludeFlag),...
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
ttlString = sprintf("Trials crossing \\theta (%s)", puffIntensity);
arrayfun(@(x) xlabel(probFigs(x).CurrentAxes, xlArray(x)), 1:Ns)
arrayfun(@(x) set(probFigs(x).CurrentAxes, axOpts{:}), 1:Ns);
arrayfun(@(x) title(probFigs(x).CurrentAxes, [ttlString;...
    "Excluded trials: "+ string(Net)+ " (\sigma:" + string(sigTh(x)) + ")"]),...
    1:Ns)
pfnStrFormat = " trial crossing PI%s TH%.1f - %.1f STH%.1f RW%.1f - %.1f s";
arrayfun(@(x) saveFigure(probFigs(x), fullfile(figureDir, sNames(x) +...
    sprintf(pfnStrFormat, puffIntensity, thSet{x}([1,end]), sigTh(x),...
    responseWindow)), 1), 1:Ns)

%% Save?
% If RollerSpeed.mat doesn't exist in the experiment folder, then create
% it.
rsfName = fullfile(dataDir, "BehaviourData" + string(expDate) + ".mat");
if ~exist(rsfName, 'file')
    save(rsfName, "genProb", "thSet", "behResults", "behStack", "stTx")
end
%% EEG-like plots for the nose and the two whisker-sets
efig = gobjects(size(behStack,1),1);
for cbp = 1:length(behStack)
    plotEEGchannels(behStack{cbp}', [], diff(timeLapse),...
        fr, 1, -timeLapse(1));
    title(sNames(cbp)); efig(cbp) = figure(); 
    plot(stTx, behStack{cbp}-median(behStack{cbp},1), "LineWidth", 0.5,...
        "Color", 0.75*ones(3,1)); hold on;
    plot(stTx, mean(behStack{cbp}-median(behStack{cbp},1),2),...
        "LineWidth", 2, "Color", "k"); title(sNames(cbp))
    lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
    figname = ['SumOfSignals', num2str(sNames(cbp)), ' filtered'];
    set(efig(cbp).CurrentAxes, axOpts{:})
    saveFigure(efig(cbp), fullfile(figureDir, figname), 1);
    
end

fprintf(1, "General probability: %.3f \n", genProb)
end
