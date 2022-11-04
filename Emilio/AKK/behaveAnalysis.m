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
axOpts = {'Box', 'off', 'Color', 'none'};
if ~exist(figureDir, "dir")
    if ~mkdir(figureDir)
        fprintf(1, "Couldn't create %s!\n", figureDir)
    end
end
%% Roller speed and trigges from both arduino and intan
% Triggers
fullname = @(x) fullfile(x.folder, x.name);
csvFiles = dir(fullfile(dataDir, 'Roller_position*.csv'));
trigFiles = dir(fullfile(dataDir, 'ArduinoTriggers*.mat'));
rsFiles = dir(fullfile(dataDir, "RollerSpeed*.mat"));
expDate = getDates(csvFiles, "Roller_position");
vars2load = {'atTimes', 'atNames', 'itTimes', 'itNames', 'Nt', 'minOfSt'};
if isempty(trigFiles)
    % Create Arduino trigger file(s)
    readAndCorrectArdTrigs(dataDir);
    trigFiles = dir(fullfile(dataDir, 'ArduinoTriggers*.mat'));
end

if isempty(rsFiles)
    % Create roller speed
    [~, vf, rollTx, fr, Texp] = createRollerSpeed(dataDir);
else
    load(fullname(rsFiles), 'vf', 'rollTx','fr','Texp')
end
% Loading the triggers
tStruct = arrayfun(@(x) load(fullname(x), vars2load{:}), trigFiles);
% Adding experiment(s) offset to arduino triggers
atT = arrayfun(@(x, z) cellfun(@(y, a) y+a, ...
    x.atTimes, repmat(z,1,length(x.atTimes)), fnOpts{:}), tStruct', ...
    num2cell([0, Texp(1:end-1)]), fnOpts{:}); atT = cat(1, atT{:});
atTimes = arrayfun(@(x) cat(1, atT{:,x}), 1:size(atT), fnOpts{:});
% Adding experiment(s) offset to intan triggers
Nt2 = [0, tStruct.Nt];
itT = arrayfun(@(x,z) cellfun(@(y, a) y+a, x.itTimes, ...
    num2cell(repmat(z,1,length(x.itTimes))), fnOpts{:}), tStruct', ...
    Nt2(1:end-1), fnOpts{:}); itT = cat(1, itT{:});
itTimes = arrayfun(@(x) cat(1, itT{:,x}), 1:size(itT), fnOpts{:});
atTimes = atTimes{1}; itTimes = itTimes{1};
% Roller speed
rollTx = rollTx{1};
% Encoder transformation to cm / s
en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fr;

%% Speed figure
fFig = figure(); 
plot(rollTx - rollTx(1), vf, "DisplayName", "Roller velocity (speed)")
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

sSig = squeeze(std(behStack(:,bsFlag,:), [], 2));
sMed = squeeze(median(behStack(:,bsFlag,:), 2));
tMed = squeeze(median(behStack, 2));

splOut = fitSpline((0:size(behDLCSignals_filtered,1)-1)/fr, behDLCSignals_filtered(:,2), 1, 0.1146, 0.9528, true);

excludeFlag = sSigma > sigTh;
% This doesn't take into consideretion the nose (1,2,4)
excludeFlag = any(excludeFlag(2.^(0:2),:),1); 
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
rsfName = sprintf("BehaviourData%s_VW%.2f - %.2f s RW%.2f - %.2f ms.mat", ...
    string(expDate), timeLapse, responseWindow*1e3);
rsfPath = fullfile(dataDir, rsfName);
if ~exist(rsfPath, 'file')
    save(rsfPath, "genProb", "thSet", "behResults", "behStack", "stTx")
end
%% EEG-like plots for the nose and the two whisker-sets
efig = gobjects(size(behStack,1),1);
yaxName = [repmat("Angle [°]",1,3), "Roller speed [cm/s]"];
for cbp = 1:length(behStack)
    plotEEGchannels(behStack{cbp}', [], diff(timeLapse),...
        fr, 1, -timeLapse(1));
    title(sNames(cbp)); efig(cbp) = figure(); 
    plot(stTx, behStack{cbp}(:,~excludeFlag)- ...
        median(behStack{cbp}(:,~excludeFlag),1), "LineWidth", 0.5,...
        "Color", 0.75*ones(3,1)); hold on;
    plot(stTx, mean(behStack{cbp}(:,~excludeFlag)- ...
        median(behStack{cbp}(:,~excludeFlag),1),2),...
        "LineWidth", 2, "Color", "k"); title(sNames(cbp))
    lgnd = legend(puffIntensity); set(lgnd, "Box", "off", "Location", "best")
    figName = sprintf("Mean %s VW%.1f - %.1f s RW%.2f - %.2f ms NEX%d", ...
        sNames(cbp), timeLapse, responseWindow*1e3, Net);
    set(efig(cbp).CurrentAxes, axOpts{:})
    xlabel(efig(cbp).CurrentAxes, "Time [s]");
    ylabel(efig(cbp).CurrentAxes, yaxName(cbp))
    saveFigure(efig(cbp), fullfile(figureDir, figName), 1);
end
fprintf(1, "General probability: %.3f \n", genProb)
close all
end
