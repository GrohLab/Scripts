data_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch18_ephys\MC\GADi43\240227_C+F_2200";
load( fullfile( data_path, "Regression CW-800.00-800.00ms " + ...
    "DW-100.00-100.00 BZ5.00.mat"), 'DX', 'mdlAll_ind', 'params')
load( fullfile( data_path, "ephys_E1", ...
    "GADi43_C+F_2200 RW20.00-200.00 SW-300.00--120.00 " + ...
    "VW-300.00-400.00 ms PuffAll (unfiltered) RelSpkTms.mat"), ...
    'relativeSpkTmsStruct' )
fnOpts = {'UniformOutput', false};
vec2trials = @(x) reshape( x, params.Nb, params.Nr, params.Ns );
bpn_abb = {'SWM', 'SWF', 'NWM', 'NWF', 'WA', 'S', 'N', 'RS'};
m = 1e-3; k = 1e3;
nox = @(x) set( get( x, 'XAxis' ), 'Visible', 'off' );
lnOpts = {'Marker', '|', 'MarkerSize', 12, 'LineStyle', 'none', 'Color', 'k'};
%% Auxiliary variables for decisions and plotting
rel_win = params.relative_window;
bin_size = params.bin_size;

mdl_mu = squeeze( mean( mdlAll_ind, 2 ) );
y_trials = vec2trials( DX{1} );
y_pred = DX{2} * mdl_mu;
y_ptrials = vec2trials( y_pred );
trial_tx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

SSEt = squeeze( sum( ( y_trials - y_ptrials ).^2, 1 ) );
SSTt = squeeze( sum( ( y_trials - mean( y_trials, 1 ) ).^2, 1 ) );
r_sq_trials = 1 - (SSEt./SSTt);

% Decide based on which trial has the highest R²: trial 6, NWM (3)
[~, ord] = sort( r_sq_trials, "descend" );
rst = relativeSpkTmsStruct(1).SpikeTimes;
%%
f = figure( "Color", "w" ); t = createtiles( f, 4, 1 );
ax = gobjects( 3, 1 );

ax(1) = nexttile( t );
line( ax(1), k*trial_tx, y_trials(:,6,3), 'LineWidth', 2, 'Color', 0.15*ones(1,3) )
cleanAxis(ax(1)); nox(ax(1));
yticklabels( ax(1), yticks(ax(1))-90 )
ylabel( ax(1), 'Angle [°]'); ytickangle( ax(1), 90 )

ax(2) = nexttile( t, [2,1] ); hold( ax(2), "on" )
Nu = size( rst, 1 );
Npu = 15;
if Npu > Nu
    Npu = Nu;
end
ct = 1;
for cu = round(linspace(1, Nu, Npu))
    if cu > Nu
        break
    end
    rst_aux = rst{cu,6};
    line( ax(2), k*rst_aux, ct+zeros( size( rst_aux, 2 ), 1 ), lnOpts{:} )
    ct = ct + 1;
end
nox(ax(2)); cleanAxis(ax(2));
ylabel(ax(2), 'Units'); ytickangle( ax(2), 90 )

ax(3) = nexttile( t );
line( ax(3), k*trial_tx, y_ptrials(:,6,3), 'LineWidth', 2, 'Color', [0,0.6,0] )
cleanAxis(ax(3));
yticklabels( ax(3), yticks(ax(3))-90 )
ylabel( ax(3), 'Angle [°]'); ytickangle( ax(3), 90 )
set( ax, 'TickDir', 'out' )
linkaxes( ax, 'x'); xlim( ax(3), [-200, 400] )
xlabel( ax(3), 'Time [ms]' )

xline( ax(3), -200:20:400, 'Color', 0.5*ones(1,3) )
