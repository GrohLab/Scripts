function [mavpt] = getMaxAbsPerTrial(inStack, responseWindow, timeAxis)
%GETMAXABSPERTRIAL gets as the name suggests, the maximum absolute
%amplitude per trial in the given input stack
%   Detailed explanation goes here, later.

%% 
% Size of the stack: Number of signals (e.g. nose, left whisker), number of
% time samples, and number of triggers (trials)
[Ns, Nts, Ntg] = size(inStack);

% Preallocating the output: Maximum Absolute Value Per Trial
mavpt = zeros(Ntg, Ns);

% Validation for time axis. If the axis doesn't correspond to the stack,
% then we cannot continue.
if Nts ~= length(timeAxis)
    fprintf(1, "Time axis not the same size as the time dimension of the stack!\n");
    fprintf(1, "Please, verify they are the same size!\n");
    return
end
% Response period for all trials
responseFlags = timeAxis >= responseWindow(1) &...
    timeAxis <= responseWindow(2);

% Loop through all signals and triggers.
for cs = 1:Ns % Signals
    for ctg = 1:Ntg % Trials or triggers
        mavpt(ctg, cs) = max(abs(inStack(cs, responseFlags, ctg)));
    end
end
end

