fnOpts = {"UniformOutput", false};
vidTx = (0:length(sPx)-1)'/fr;
%expTx = (0:length(trig)-1)'/fs;

timeLapse = [-2, 4];
consideredConditions = 3:5;
% All triggers in the experiment without overlap
allTSubs = sort(cat(1,Conditions(consideredConditions).Triggers));
% Aligning the roller signal to the highest temporal correlation.
[xcr, lgs] = xcorr(sPx, vf); [~, mxLg] = max(xcr); delSamples = lgs(mxLg);
del = delSamples / fsRoll;
vf_corrected = cat(1, zeros(delSamples, 1), vf(:));
% Creating the stack for the roller speed. Assumption: Acustic stimulation
% makes the mouse to walk backwards (or move the roller in any direction).
[~, vStack] = getStacks(false, allTSubs, 'on', timeLapse, fs, fsRoll,...
    [], {vf_corrected}); [Ne, Nt, NTa] = size(vStack);
vStack = squeeze(vStack);
[ms, bs] = lineariz([1,size(vStack,1)], timeLapse(2), timeLapse(1));
% Stack time axis
stTx = (1:size(vStack,1))' * ms + bs; spontFlag = stTx < 0;
excludeFlag = rms(vStack(spontFlag,:)) > 0.85 | rms(vStack) > 0.9;
plotEEGchannels(vStack(:,~excludeFlag)', [], diff(timeLapse), fsRoll,...
    1, abs(timeLapse(1)));
% Peaks and valleys for each trial.
cpts = getWaveformCriticalPoints(vStack, fsRoll);
cpts = cellfun(@(x) x + timeLapse(1), cpts, fnOpts{:});
% Amplitude for each peak and valley
camp = arrayfun(@(x) interp1(stTx, vStack(:,x), cpts{x,1}),...
    (1:size(vStack,2))', fnOpts{:});
% Location for the critical point
[~, mxSub] = cellfun(@max, camp, fnOpts{:});
dlTms = arrayfun(@(x) cpts{x,1}(mxSub{x}), (1:size(cpts,1))', fnOpts{:});
ttms = allTSubs./fs; emptyCells = cellfun(@isempty, dlTms);
% For a given experimental time i.e. trigger time we estimate an assosiated
% delay for the Arduino board.
pts = [ttms(~emptyCells), cat(1,dlTms{~emptyCells})];
figure; scatter(pts(:,1), pts(:,2), "filled")
%% DE_jittering
Nccond = length(consideredConditions); 
delayFlags = false(NTa,Nccond); chCond = 1;
counter2 = 1;
for ccond = consideredConditions
    delayFlags(:,counter2) = ismember(allTSubs(:,1),...
        Conditions(ccond).Triggers(:,1)); counter2 = counter2 + 1;
end
Na = sum(delayFlags,1);
%% RANSAC line estimation
n = 1;
[rmdl, inln] = boot_fit_poly([pts(:,1), pts(:,2)-0.5], n, (n+1)/size(pts,1),...
    2^17, 0.56);
%% RANSAC correction
Conditions_corrected = arrayfun(@(x) struct('name', Conditions(x).name,...
    'Triggers', Conditions(x).Triggers +...
    round((Conditions(x).Triggers./fs) .* rmdl(1) + rmdl(2))*fs),...
    (1:length(Conditions))');

%% RANSAC threshold estimation
n = 1; xpts = [0;pts(end,1)];
for cth = 0.5:0.01:0.6
    [rmdl, inln] = boot_fit_poly(pts, n, (n+1)/size(pts,1), 2^17, cth);
    if islogical(inln)
        figure; gscatter(pts(:,1), pts(:,2), inln); title(sprintf('%.4f',cth))
        hold on; plot(xpts, xpts.^(n:-1:0) * rmdl)
    end
end
%% Special case for GAD49 frequency stimulation
Conditions_corrected = arrayfun(@(x) struct('name', Conditions(x).name, 'Triggers', Conditions(x).Triggers + ...
    round(cat(1, (Conditions(x).Triggers(Conditions(x).Triggers(:,1)<(675.8*fs),:)/fs) .* mdl_1(1) + mdl_1(2) - 1, ...
    (Conditions(x).Triggers(Conditions(x).Triggers(:,1) >= (675.8*fs) & Conditions(x).Triggers(:,1) < (887.3*fs),:)/fs) .* mdl_2(1) + mdl_2(2) - 1.3, ...
    (Conditions(x).Triggers(Conditions(x).Triggers(:,1) >= (887.3*fs) & Conditions(x).Triggers(:,1) < size(expTx,1),:)/fs) .* mdl_3(1) + mdl_3(2) - 1)*fs)),...
    (1:length(Conditions))');

%% Producing the final figures
en2cm = ((2*pi)/((2^15)-1))*((14.85/2)^2)*fsRoll;
miinOpts = {"Color", 'k', "LineWidth", 2};
vMean = zeros(length(stTx), Nccond); vcount = 1;
vFigs = gobjects(Nccond+1, 1); 
for ccond = consideredConditions
    [~, vStack] = getStacks(false, Conditions_corrected(ccond).Triggers, 'on',...
        timeLapse, fs, fsRoll, [], {vf_corrected}); vStack = squeeze(vStack);
    % No movement before nor after
    excludeFlag = rms(vStack(spontFlag,:)) > 0.85 | rms(vStack) > 0.9;
    % Movement before
    % excludeFlag = rms(vStack(spontFlag,:)) < 0.85;
    % No exclusion
    % excludeFlag = false(size(vStack,2),1);
    % No movement before but only after the trigger
    % excludeFlag = rms(vStack(spontFlag,:)) > 0.85 | rms(vStack) < 0.9;
    % No movement before
    % excludeFlag = rms(vStack(spontFlag,:)) > 0.85;
    if all(excludeFlag)
        continue
    end
    plotEEGchannels(vStack(:, ~excludeFlag)', '', diff(timeLapse),...
        fsRoll, 1, abs(timeLapse(1)));
    vFigs(vcount) = figure; 
    plot(stTx, vStack(:, ~excludeFlag)*en2cm, "Color", 0.65*ones(3,1));
    vMean(:,vcount) = mean(vStack(:, ~excludeFlag),2)*en2cm;
    hold on; plot(stTx, vMean(:,vcount), miinOpts{:})
    set(gca, "Box", "off", "Color", "none", "Clipping", "off"); 
    title(Conditions(ccond).name)
    xlabel("Time [s]"); ylabel("Roller speed [cm/s]")
    vcount = vcount + 1;
end
%%
arrayfun(@(x) saveFigure(vFigs(x), fullfile(figureDir,...
    sprintf('RollerSpeed_%s VW%.1f - %.1f s (no movement)',...
    Conditions(x+2).name, timeLapse)),1, 1), (1:3)');
vFigs(4) = figure; plot(stTx, vMean); 
lgnd = legend(arrayfun(@(x) Conditions(x).name, consideredConditions, fnOpts{:}));
set(lgnd, "Box", "off", "Location", "best")
set(gca, "Box", "off", "Color", "none")
xlabel("Time [s]"); ylabel("Roller speed [cm/s]")
title("Condition comparison")
saveFigure(vFigs(4), fullfile(figureDir,...
    sprintf(...
    'Roller mean speed comparison_3 conditions VW%.1f - %.1f s (no movement)',...
    timeLapse)),1, 1);
%%
figure;
for cm = unique(mSubs)'
    scatter3(coeff(mSubs == cm, 1), coeff(mSubs == cm, 2),...
        coeff(mSubs == cm, 3), "filled", "DisplayName", num2str(cm),...
        "MarkerFaceAlpha", 0.8)
    if cm == 1
        hold on
    end
end