function [outputArg1,outputArg2] = AnBeh_Bypass(data_path, resWin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% Input validation

checkDir = @(x) (ischar(x) || isstring(x)) && exist(x, 'dir');
checkRW = @(x) isnumeric(x) & numel(x) == 2 & (x(2) > x(1));

p = inputParser;

addRequired(p, 'data_path', checkDir)
addRequired(p, 'resWin', checkRW)

parse(p, data_path, resWin);

data_path = p.Results.data_path;
brWin = p.Results.resWin;
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};
createtiles = @(f,r,c) tiledlayout( f, r, c, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
lgOpts = {'Box', 'off', 'Color', 'none', 'Location', 'best'};
cleanAxis = @(x) set( x, lgOpts{1:4} );
vaxOpts = cellstr( ["HorizontalAlignment", "center", ...
    "VerticalAlignment", "baseline", "Rotation"] );
%TODO: Get bodypart names from somewhere else!
bodypart_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
"Nonstim-whisker mean", "Nonstim-whisker fan arc", "Interwhisk arc", ...
"Symmetry", "Nose", "Roller speed"];
expandName = @(x) fullfile( x.folder, x.name );
m = 1e-3; k = 1e3;
%% Load regression file
regFile = dir( fullfile( data_path, "Regression CW*.mat" ) );
if ~isempty( regFile ) && numel( regFile ) == 1
    load( expandName( regFile ), 'DX', 'mdlAll_ind', 'params')
else
    %TODO: Being able to select a regression for specific parameters
    fprintf(1, 'Regression either didn''t work or have different files\n')
    return
end
% Colormap: grey for original and PCB-green for reconstructed
clrMap = flip([0.15*ones(1,3); 0, 51/255, 0], 1);
% Getting parameter values
rel_win = params.relative_window; % [-1, 1]*0.8;
bin_size = params.bin_size; % 10*m;
Nb = params.Ns;
% Computing time axis for trials
trial_tx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

% Spontaneous window
bsWin = -flip(brWin);

brWin_aux = brWin;
if brWin(1) < 0.12
    brWin_aux = brWin + 0.1;
end
brWin = [brWin_aux; repmat( brWin, Nb-1, 1 )];


%% Organising figures in subfolders
vwKey = sprintf("V%.2f - %.2f s", rel_win);
dwKey = sprintf("D%.2f - %.2f ms", del_win*k);
rwKey = sprintf("R%.2f - %.2f ms", brWin(end,:)*k);
bsKey = sprintf("B%.2f ms", bin_size*k);
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
vec2tr = @(x) reshape( x, params.Nb, [], params.Ns );

% Reconstructing behaviour (Laser OFF only)
mdl_mu = squeeze( mean( mdlAll_ind, 2 ) );
rb = cellfun(@(x) x * mdl_mu, DX([2,3]), fnOpts{:} );
% Reorganising to have trials
rbStack = cellfun(vec2tr, rb, fnOpts{:} );
obStack = cellfun(vec2tr, DX([1,4]), fnOpts{:} );
[Nts, Nr, Nb] = size( rbStack{1} );
stk = cat( 1, rbStack(1), obStack(1) );
mvpt = zeros( Nb, Nr, 2 );
up_pc = 1.15;
for cs = 1:2
    % Getting maximum speed per trial
    mvpt_aux = arrayfun(@(b) getMaxAbsPerTrial( squeeze( stk{cs}(:,:,b) ), ...
        brWin(b,:), bsWin, trial_tx ), 1:Nb, fnOpts{:} );
    mvpt(:,:,cs) = cat( 1, mvpt_aux{:} );
end
% Normalising each
mvps = max( mvpt, [], 2 ) * up_pc;
ai = squeeze( mean( mvpt ./ mvps, 2 ) );

radAxis = (0:Nb-1)*(2*pi/Nb);
z_axis = exp(1i*radAxis(:));

poly_coords = ai .* z_axis;

f = figure('Color', 'w'); t = createtiles( f, 1, 1 ); ax = nexttile(t);
% Drawing polar axis 
line(ax, [-real(z_axis(1:4)), real(z_axis(1:4))]', ...
    [-imag(z_axis(1:4)), imag(z_axis(1:4))]', 'LineWidth', 0.1, ...
    'Color', 0.45*ones(1,3));
arrayfun(@(x) rectangle(ax, 'Position', [repmat(-x,1,2), repmat(x*2,1,2)], ...
    'Curvature', [1,1], 'EdgeColor', 0.45*ones( 1, 3 ) ), 1:-0.25:0.25 );
text( ax, 0.25:0.25:1, zeros(1,4), string( ( 1:4 )'/4 ), ...
    "HorizontalAlignment", "left", "VerticalAlignment", "cap" )

pchObj = arrayfun(@(c) patch( ax, real( poly_coords(:,c) )', ...
    imag( poly_coords(:,c) )', clrMap(c,:), 'EdgeColor', 'none', ...
    'FaceAlpha', 0.5 ), 1:2 );

legend( ax, pchObj, {'Reconstructed', 'Observed'}, lgOpts{:} );

arrayfun(@(v,b,y) text( ax, real( z_axis(v) ), ...
    imag( z_axis(v) ), b, vaxOpts{:}, y ), 1:Nb, bodypart_names, ...
    (180*angle( transp( z_axis ) )/pi) - 90 )

axis(ax, 'off', 'equal'); cleanAxis( ax );



end