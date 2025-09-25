function [mavpt, mxT, bl_lvl] = getMaxAbsPerTrial(inStack, ...
    responseWindow, spontWindow, timeAxis)
%GETMAXABSPERTRIAL gets as the name suggests, the maximum absolute
%amplitude per trial in the given input stack
%   Detailed explanation goes here, later.

%% 
my_xor = @(x) xor( x(:,1), x(:,2) );
% Size of the stack: Number of signals (e.g. nose, left whisker), number of
% time samples, and number of triggers (trials)
[Nts, Ntg] = size(inStack);

% Preallocating the output: Maximum Absolute Value Per Trial
mavpt = zeros(Ntg, 1);

% Validation for time axis. If the axis doesn't correspond to the stack,
% then we cannot continue.
if Nts ~= length(timeAxis)
    fprintf(1, "Time axis not the same size as the time dimension of the stack!\n");
    fprintf(1, "Please, verify they are the same size!\n");
    return
end
% Response and spontaneous periods for all trials
responseFlags = my_xor( timeAxis(:) >= responseWindow );
spontaneousFlags = my_xor( timeAxis(:) >= spontWindow );

% mavpt = max(abs(inStack(responseFlags, :)-median(inStack,1)));
bl_lvl = median( inStack(spontaneousFlags,:), 1 );
[mavpt, ps] = max( abs( inStack(responseFlags,:) - bl_lvl ) );
% mavpt = arrayfun(@(mp, tr) ...
%     inStack(find(responseFlags,1,'first')+mp-1, tr), ps(:), (1:Ntg)');
mxT = timeAxis(ps+find(responseFlags,1,"first")-1);
end