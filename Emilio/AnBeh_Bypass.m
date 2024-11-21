function [outputArg1,outputArg2] = AnBeh_Bypass(inputArg1,inputArg2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rel_win = params.relative_window; % [-1, 1]*0.8;
del_win = params.delay_window; % [-100, 100]*m;
bin_size = params.bin_size; % 10*m;

% Spontaneous flags
bsWin = -flip(brWin);
bsFlag = behTx < bsWin; bsFlag = xor(bsFlag(:,1), bsFlag(:,2));

brWin_aux = brWin;
if brWin(1) < 0.12
    brWin_aux = brWin + 0.1;
end

brWin = [brWin_aux; repmat( brWin, Nbs-1, 1 )];

brFlag = arrayfun(@(x) behTx < brWin(x,:), 1:Nbs, fnOpts{:} );
brFlag = cellfun(@(x) xor(x(:,1), x(:,2) ), brFlag, fnOpts{:} );
brFlag = cat( 2, brFlag{:} );

%% Organising figures in subfolders
vwKey = sprintf("V%.2f - %.2f s", bvWin);
rwKey = sprintf("R%.2f - %.2f ms", brWin(2,:)*k);
subFig = "Beh %s %s";
% Configuration subfolder
subfigDir = fullfile(figureDir, sprintf(subFig, vwKey, rwKey));
metaNameFlag = false;
if exist(subfigDir, "dir")
    figureDir = subfigDir;
else
    if ~mkdir(subfigDir)
        % Print metadata on figure name
        metaNameFlag = true;
        if verbose
            fprintf(1, "Error while creating subfolder!\n")
            fprintf(1, "Placing the figures in 'Figure' directory.\n")
        end
    else
        figureDir = subfigDir;
    end
end

trial_tx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);


% Getting maximum speed per trial
mvpt = arrayfun(@(c) arrayfun( @(b) ...
    getMaxAbsPerTrial( behStack{b}(:, xtf(:, b, c)), brWin(b,:), behTx ), ...
    1:Nbs, fnOpts{:} ), 1:Nccond, fnOpts{:} );
mvpt = cat(1, mvpt{:});
end