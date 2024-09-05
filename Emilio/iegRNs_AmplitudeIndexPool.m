
fnOpts = {'UniformOutput', false};
exp_subtype = {'iRNs', 'eRNs', 'RNs'};
% exp_subtype = {'terminal inhib'};
% name_keys = {'GADi', {'GADe', 'vGlut'}, 'WTg'};
% name_keys = {'GADi', 'vGlut', 'WTg'};
name_keys = {'GADi', 'GADe', 'WTg'};
% name_keys = {{'Rb', 'WT'}};
bodypart_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", "Interwhisk arc", ...
    "Symmetry", "Nose", "Roller speed"];
signTh = [0.001, 0.01, 0.05, 0.1];
createtiles = @(f,nr,nc) tiledlayout( f, nr, nc, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
bxOpts = cellstr(["Notch", "on", "JitterOutliers", "on", ...
    "BoxFaceColor", "k", 'MarkerStyle', '.', 'MarkerColor', 'k']);
txOpts = {'HorizontalAlignment','left','VerticalAlignment', 'middle', ...
    'Rotation', 90};
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );
fig_path = fullfile( ...
    "Z:\Emilio\SuperiorColliculusExperiments\Roller\PoolFigures" );
load( fullfile( fig_path, "MC, BC, BS, MCterminals pool.mat" ), "summMice" );

expMice = summMice{1}(3);
ovwtFlag = false;
% xLabels = ["Control", "C100", "F100"];
% xLabels = ["Control", "C30", "F30", "C100", "F100", ...
% "C400", "F400", "C600", "Musc", "Musc" ];
xLabels = ["Control", "C30", "F30", "C100", "F100", ...
    "C400", "F400", "Dead", "PTX", "PTX" ];
exp_subtype_flags = cellfun(@(x) contains( expMice.MiceNames, x ), ...
    name_keys, fnOpts{:} );
exclude_names = {'GADi13', 'GADi15', 'GADi53'};
exclude_mice = contains(expMice.MiceNames, exclude_names);
exp_subtype_flags = cat( 2, exp_subtype_flags{:} );
figs = gobjects( numel( exp_subtype ), 1 );
% exp_type = join( ['MC-', expMice.ExperimentalGroup] );
exp_type = expMice.ExperimentalGroup;
for cest = 1:numel(exp_subtype)
    cons_mice = exp_subtype_flags(:,cest) & ~exclude_mice;
    exp_subtype_flags(:,cest) = exp_subtype_flags(:,cest) & ~exclude_mice;
    if sum( cons_mice )

        figs(cest) = figure( "Color", "w" );
        t = createtiles( figs(cest), 2, 1 ); ax = nexttile(t);
        aux = squeeze( mean( expMice.AmplitudeIndex(:,:, ...
            cons_mice), 1, "omitmissing" ) )';
        boxchart( ax, aux, bxOpts{:} );
        ylabel( ax, 'Amplitude index' ); xticklabels( ax, xLabels );
        ylim(ax, [0,1]); cleanAxis(ax); ax.XAxis.Visible = "off";
        p = arrayfun(@(x) signrank( aux(:,1), aux(:,x) ), 2:size(aux,2), ...
            "ErrorHandler", @(s,a) nan(1) );
        mark_flag = p(:) < signTh;
        text( ax, 2:size(aux,2), 1.15*max( aux(:,2:end), [], 1 ), ...
            arrayfun(@(x) sprintf("p=%.3f", x), p(:) ), txOpts{:} )

        ax = nexttile(t);
        aux = squeeze( mean( expMice.TrialProportions(:,:, ...
            cons_mice), 1, "omitmissing" ) )';
        boxchart(ax, aux, bxOpts{:} );
        ylabel( ax, 'Trial proportions' )
        xticklabels( xLabels ); ylim([0,2]); cleanAxis( ax );
        p = arrayfun(@(x) signrank( aux(:,1), aux(:,x) ), 2:size(aux,2), ...
            "ErrorHandler", @(s,a) nan(1) );
        text(ax, 2:size(aux,2), 1.15*max( aux(:,2:end), [], 1 ), ...
            arrayfun(@(x) sprintf("p=%.3f", x), p(:) ), txOpts{:} )
        title(t, sprintf( "Area/%s", exp_subtype{cest} ) )

        for cbp = 1:numel(bodypart_names)
            fig = figure("Color", "w"); t2 = createtiles( fig, 2, 1 );
            ax = nexttile(t2);
            aux = squeeze( mean( expMice.PolygonUnfoldAmplIndx(:,:, ...
                cbp, cons_mice), 2, "omitmissing" ) )';
            boxchart( ax, aux, bxOpts{:} );
            ylabel( ax, 'Amplitude index' ); xticklabels( ax, xLabels );
            ylim(ax, [0,1]); cleanAxis(ax); ax.XAxis.Visible = "off";
            p = arrayfun(@(x) signrank( aux(:,1), aux(:,x) ), 2:size(aux,2), ...
                "ErrorHandler", @(s,a) nan(1) );
            text( 2:size(aux,2), 1.15*max( aux(:,2:end), [], 1 ), ...
                arrayfun(@(x) sprintf("p=%.3f", x), p(:) ), txOpts{:} )

            ax = nexttile(t2);
            aux = squeeze( mean( expMice.PolygonUnfoldTrialProp(:,:, ...
                cbp, cons_mice), 2, "omitmissing" ) )';
            boxchart(ax, aux, bxOpts{:} );
            ylabel( ax, 'Trial proportions' )
            xticklabels( xLabels ); ylim([0,2]); cleanAxis( ax );
            p = arrayfun(@(x) signrank( aux(:,1), aux(:,x) ), 2:size(aux,2), ...
                "ErrorHandler", @(s,a) nan(1) );
            text( 2:size(aux,2), 1.15*max( aux(:,2:end), [], 1 ), ...
                arrayfun(@(x) sprintf("p=%.3f", x), p(:) ), txOpts{:} )
            title(t2, sprintf( "%s/%s", bodypart_names(cbp), exp_subtype{cest} ) )

            % saveFigure( fig, fullfile( fig_path, join( [bodypart_names(cbp), ...
            %     exp_type, exp_subtype{cest}, "all mice pool" ] ) ), ...
            %     true, ovwtFlag )
            saveFigure( fig, fullfile( fig_path, join( [bodypart_names(cbp), ...
                exp_type, exp_subtype{cest}, sum( cons_mice ) ] ) ), ...
                true, ovwtFlag )
        end
    end
end

% arrayfun(@(f) saveFigure( figs(f), fullfile( fig_path, ...
%     join( ["Areas", exp_type, exp_subtype{f}, "all mice pool"] ) ), true, ovwtFlag ), ...
%     find( arrayfun(@(f) ~isa( f, 'matlab.graphics.GraphicsPlaceholder'), figs ) ) );

arrayfun(@(f) saveFigure( figs(f), fullfile( fig_path, ...
    join( ["Areas", exp_type, exp_subtype{f}, ...
    sum( exp_subtype_flags(:,f) ) ] ) ), true, ovwtFlag ), ...
    find( arrayfun(@(f) ~isa( f, 'matlab.graphics.GraphicsPlaceholder'), figs ) ) );