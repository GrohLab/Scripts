%%
mice_results = "Z:\Emilio\SuperiorColliculusExperiments\Roller\PoolFigures\MC-iegRNs";
pool_fig_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller\PoolFigures\MC-iegRNs\iRNs";

load( fullfile( mice_results, 'iRNs\MCiRNs_reconstruction_sm.mat') )
bp_names = ["Stim-whisker mean", "Stim-whisker fan arc", ...
    "Nonstim-whisker mean", "Nonstim-whisker fan arc", ...
    "Interwhisker arc", "Symmetry", "Nose", "Roller speed"];

Nm = numel( mice ); % Numer of mice
Nspm = arrayfun(@(x) numel( x.Sessions ), mice ); % Number of sessions per mouse
Nexp = sum( Nspm );
r2_res_c = zeros( 3, 8, Nexp );
r2_res_l = zeros( 3, 8, Nexp );
ce = 1;
mouseID = zeros( Nexp, 1 );
sessID = zeros( Nexp, 1 );
for cm = 1:Nm
    for cs = 1:Nspm(cm)
        dt = mice(cm).Sessions(cs).DataTable;
        r2 = dt.R_2_p_L;
        r2_res_c(:,:,ce) = r2{1};
        r2_res_l(:,:,ce) = r2{2};
        mouseID(ce) = cm;
        sessID(ce) = cs;
        ce = ce + 1;
    end
end
[Nep, Ns] = size( r2_res_c, [1,2] );

%%
fnOpts = {'UniformOutput', false};
tocol = @(x) x(:);
r2_mean_c = cellcat( arrayfun(@(x) mean( r2_res_c(:,:,mouseID==x), 3 ), ...
    1:Nm, fnOpts{:} ), 3 );
r2_mean_l = cellcat( arrayfun(@(x) mean( r2_res_l(:,:,mouseID==x), 3 ), ...
    1:Nm, fnOpts{:} ), 3 );
f = figure("Color", "w"); t = createtiles( f, 1, 1); 
ax = nexttile( t );
bpID = repmat( ones( Nep, 1 ) * (1:Ns), 1, 1, Nm);
preVSpostID = repmat( (1:Nep)' * ones( 1, Ns ), 1, 1, Nm );
boxchart( ax, bpID(:), tocol( r2_mean_c ), 'GroupByColor', preVSpostID(:), ...
    'Notch', 'on' )
xline( ax, (1:Ns-1) + 1/2, '--', 'Color', 0.45*ones(1,3) ); 
legend( {'Overall', 'Pre', 'Post'}, "Box", "off", "Color", "none", ...
    "Location", "best", "AutoUpdate", "off" )
cleanAxis( ax ); ytickangle( ax, 90 ); set( ax, 'TickDir', 'out' );
ylabel( ax, 'R²' )
xticks( ax, 1:Ns ); xticklabels( ax, bp_names ); 
xlim( ax, [1,Ns] + [-1,1]/2 ); 
ylim( ax, [0, 1] )
set( f, 'UserData', {r2_mean_c, bpID, preVSpostID, mouseID, sessID} )
title( ax, 'Reconstruction R² for overall, pre-, and post-stimulus' )

%%
p = arrayfun(@(x) signrank( squeeze( r2_mean_c(2,x,:) ), ...
    squeeze( r2_mean_c(3,x,:) ) ), 1:Ns );
fnOpts = {'UniformOutput', false};
txOpts = {'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'};
astk = sum( p < [0.05, 0.01, 0.001]' );
x = (1:Ns) + [0;1]/3;
y = [1;1] * max( r2_mean_c(2:3,:,:), [], [3,1] ) * 1.05;
line( ax, x, y, 'Color', 'k' )
txt = [arrayfun( @(a) replace( join( repmat("\ast", 1, a) ), " ", "" ), astk, fnOpts{:} );
arrayfun(@(h) sprintf( "$p=%.3f$", h) , p)];
%txt = arrayfun(@(s) replace( join( txt(:,s) ), " ", ""), 1:Ns );
text( ax, mean( x, 1 ), y(1,:)+0.035, txt(1,:), txOpts{:}, "FontSize", 10 )
text( ax, mean( x, 1 ), y(1,:), txt(2,:), txOpts{:}, "FontSize", 8, ...
    "Interpreter", "latex" )
%%
saveFigure( f, fullfile( pool_fig_path, ...
    "Reconstruction R² overall, pre and post" ), true, true )