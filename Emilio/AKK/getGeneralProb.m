function [genProb, probMat] = getGeneralProb(stack, theta, en2cm, rmsTh1, rmsTh2)
%GETGENERALPROB takes in the stack and the threshold and returns the
%general probability of the relevant stack and its probability matrix.
% en2cm is just for the speed.

%The Problem here is that I do not know how I distinguish between different stacks
%as input arguments.

% cell- and arrayfun auxiliary variable.
fnOpts = {'UniformOutput', false};

%Not sure if these are the same or if they have to be fitted:
timeLapse = [-1, 2];
responseWindow = [0, 0.4]; % will stay the same right?

%% if loop: if stack is speed or movement ??? or just one for all? 
%I could also do one function for speedprob, one for the bodyparts etc.
%Right now I write it like it's just one for all:

[Nts, NTa] = size(stack);
[ms, bs] = lineariz([1, Nts], timeLapse(2), timeLapse(1));

% Stack time axis --> st T x
stTx = (1:Nts)'*ms + bs;
% Spontaneous logical flag.
spontFlag = stTx < 0;
% No previous movement before nor after the stimulus
% Logical Flag to exclude trials with 'too much movement' before or after
% the trigger.
excludeFlag = rms(stack(spontFlag,:)) > rmsTh1 | rms(stack) > rmsTh2;

delayFlags = true(NTa,1); 
Nccond = 1;
% Number of triggers (Number of alignment points) Na excluding the
% indicated trials.
Na = sum(delayFlags & ~excludeFlag(:));

% Loop to get the maximum value of the movement.
cellMax = cell(Nccond,1);
responseFlags = stTx >= responseWindow(1) & stTx <= responseWindow(2);

for ccond = 1:Nccond
    cellMax{ccond} =...
        max(abs(stack(responseFlags, delayFlags(:, ccond) & ~excludeFlag(:))))';
end

% Maximum number of trials
maxNt = max(cellfun(@numel, cellMax));
% Placing the trial's maximum movement (of any kind) in a matrix
moveMat = cellfun(@(x) cat(1, x, nan(maxNt - size(x,1),1)),...
    cellMax, fnOpts{:}); 
moveMat = cat(2, moveMat{:});

%just for speed: ?
moveMat = moveMat * en2cm;

% Cutting on different thresholds for all trials
moveFlag = moveMat > theta;
probMat = sum(moveFlag)./Na; 
probMat = squeeze(probMat);

genProb = nnz(moveFlag)./numel(moveFlag);

disp('General Probability is calculated.');

end