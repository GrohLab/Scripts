mID = 1;
Nm = numel( reg_mice );
Ns = sum( arrayfun(@(m) numel( m.Sessions ), reg_mice ) );
rst_reg_dt = zeros( Ns, 5 );
rsq_dt = zeros( Ns, 10 );
regMiceNames = [reg_mice.Name]';
rstMiceNames = [rst_mice.Name]';
cr = 1;
for cm = 1:numel(reg_mice)
    sID = 1;
    for cs = 1:numel(reg_mice(cm).Sessions)
        dt = rst_mice(cm).Sessions(cs).DataTable;
        rst_reg_dt(cr,:) = [cm, cs, dt.NUnits, dt.Proportion];
        dt = reg_mice(cm).Sessions(cs).DataTable;
        rsq_dt(cr,:) = [cm, cs, dt.R_squared];
        cr = cr + 1;
    end
end

%%

f = figure( "Color", "w" );
createtiles = @(f,r,c) tiledlayout( f, r, c, 'TileSpacing', 'Compact', ...
    'Padding', 'tight');
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );
t = createtiles( f, 3, 1 );

ax = gobjects(3, 1); pl = ax;
clrMap = [0,0,0; 0.45*ones(1,3); 0,0,0];
x = 0:0.05:0.8;
samples = [1e4, 1];
ylbls = ["Recorded units", "Responsive 20 - 50 ms", "Responsive 50 - 200 ms"];
for cp = 1:3
    ax(cp) = nexttile( t );
    pl(cp) = plot( ax(cp), rsq_dt(:,3), rst_reg_dt(:,2+cp), 'Marker', '.', ...
        'Color', clrMap(cp,:), 'MarkerSize', 16, 'LineStyle', 'none' );
    cleanAxis( ax(cp) );
    ylabel( ax(cp), ylbls(cp) );
    ft = fitlm( rsq_dt(:,3), rst_reg_dt(:, 2+cp ), 'poly1' );
    cCI = coefCI( ft ); sig = (ft.Coefficients{:,"Estimate"} - cCI(:,1))/2;
    slope = makedist( "Normal", "mu", ft.Coefficients{"x1","Estimate"}, ...
        "sigma", sig(2) );
    intercept = makedist( "Normal", "mu", ...
        ft.Coefficients{"(Intercept)","Estimate"}, "sigma", sig(1) );
    y = random( slope, samples ) * x + random( intercept, samples );
    CI = quantile( y, [0.025, 0.975], 1 );
    hold( ax(cp), 'on' ); lobjs = line( x, CI', 'Color', 0.35*ones(1,3), ...
        'LineStyle', '--' ); 
    lobjs(1) = [];
    lobjs = cat(1, lobjs, line( x, predict( ft, x' ), 'Color', 0.35*ones(1,3) ));
    legend( lobjs, {'CI 95%', sprintf('Fit %.3f r^2', ...
        ft.Rsquared.Ordinary)}, Box="off", Color='none', Location='best' );
end
xlabel( ax(cp), "R^2_{SWM}" )
set( ax, 'TickDir', 'out' )
arrayfun(@(x) set( get( x, 'XAxis' ), 'Visible', 'off' ), ax(1:2) )
%%
Bsel = 1;
Nex = 3;
f = figure("Color", "w" );
t = createtiles( f, Nex+2, 1 );

getSEM = @(x) [mean( x, 2 ), std( x, 1, 2 )./sqrt( size( x, 2 ) )];
mat2ptch = @(x) [x(1:end,:)*[1;1]; x(end:-1:1,:)*[1;-1]];
phOpts = {'EdgeColor', 'none', 'FaceAlpha', 0.25, 'FaceColor'};

clrMap = brighten( blues(2), -0.8 );
ax = gobjects( Nex+1, 1 );
[~, ord] = sort( rmse_ltrials, "ascend" );
cpi = 1;
for cp = ord(1:round(end/Nex):end)
    ax(cpi) = nexttile( t ); 
    set( ax(cpi), 'NextPlot', 'add' ); cleanAxis( ax(cpi) );
    line( ax(cpi), tr_tx, squeeze( y_ltrials(:,cp,Bsel) ), 'Color', 0.15*ones(1,3), 'LineWidth', 0.5 )
    line( ax(cpi), tr_tx, squeeze( y_lptrials(:,cp,Bsel) ), 'Color', 0.55*ones(1,3), 'LineWidth', 0.5 )
    yticklabels( ax(cpi), yticks(ax(cpi)) - 90)
    ax(cpi).YTickLabelRotation = 90;
    set( get( ax(cpi), 'XAxis' ), 'Visible', 'off' )
    cpi = cpi + 1;
end
ax(cpi) = nexttile( t, [2, 1] );
cleanAxis( ax(cpi) )
sem_lshadow = getSEM( y_ltrials );
sem_lpshadow = getSEM( y_lptrials );
line(ax(cpi), tr_tx, squeeze( mean( y_ltrials(:,:,Bsel), 2 ) ), 'LineWidth', 2, 'Color', clrMap(1,:) )
line(ax(cpi), tr_tx, squeeze( mean( y_lptrials(:,:,Bsel), 2 ) ), 'LineWidth', 2, 'Color', clrMap(2,:) )
patch( ax(cpi), [tr_tx(:); flip( tr_tx(:) )], mat2ptch( sem_lshadow(:,:,Bsel) ), 1, phOpts{:}, clrMap(1,:) )
patch( ax(cpi), [tr_tx(:); flip( tr_tx(:) )], mat2ptch( sem_lpshadow(:,:,Bsel) ), 1, phOpts{:}, clrMap(2,:) )
yticklabels( ax(cpi), yticks(ax(cpi)) - 90); ylabel( ax(cpi), 'Angle [°]')
xticklabels( ax(cpi), xticks( ax(cpi) )*1e3 );
xlabel(ax(cpi), 'Time [ms]')
legend( ax(cpi), flip({'Predicted', 'Observed'}), 'Box', 'off', 'Color', 'none', ...
    'Location', 'best', 'AutoUpdate', 'off' );
set( ax, 'TickDir', 'out' )
xline(ax(cpi), 0, 'k--')
xline( ax(cpi), [-0.1, 0.2], 'b')

saveFigure( f, fullfile( "Z:\Emilio\SuperiorColliculusExperiments\" + ...
    "Roller\Batch7_ephys\MC\GADi52\220808_C+F_2100\ephys_E1\Figures", ...
    "Example reconstruction laser" ), true, false )

%%
f = figure("Color", "w");
t = createtiles( f, 3, 1);
my_xor = @(x) xor( x(:,1), x(:,2) );
lnOpts = {'LineStyle', 'none', 'Marker', '|', 'Color', 0.15*ones(1,3)};
tlSelect = 9;
vWin = [-0.3,0.5];
[~, ord2] = sort( cellfun(@(x) size(x, 1), spike_times ), "descend" );
cni = 1;

ax = gobjects( 2, 1 );
ax(1) = nexttile( t, [2, 1] ); cleanAxis( ax(1) ); ylabel( ax(1), 'Units')
for cn = ord2(1:round(end/30):end)'
    spkIdx = my_xor( spike_times{cn} > time_limits(tlSelect,:) );
    line( spike_times{cn}(spkIdx), cni+zeros(sum(spkIdx),1), lnOpts{:} )
    cni = cni + 1;
end
set( get( ax(1), 'XAxis' ), 'Visible', 'off' )
bin_size = 0.02;
ax(2) = nexttile(t);
behIdx = my_xor( btx > time_limits(tlSelect,:) );
line(ax(2), btx(behIdx), behSignals( behIdx, 1 ) )
linkaxes( ax, 'x' )
xline(ax(2), time_limits(tlSelect,1)+(bin_size/2):bin_size: ...
    time_limits(tlSelect,2)-(bin_size/2), 'LineWidth', 0.1, 'Color', 0.5*ones(1,3) )
xlim( ax(2), mean(time_limits(tlSelect,:)) + vWin)
xline( ax(2), mean( time_limits(tlSelect,:)), 'k--' )
set( ax, 'TickDir', 'out' )
xtv = mean( time_limits(tlSelect,:) ) + vWin;
xticks( ax(2), xtv(1):0.1:xtv(2) ); xlabel( ax(2), 'Time [ms]' )
xticklabels( ax(2), round( ( xticks( ax(2) ) - xtv(1) + vWin(1) ) *1e3 ) )
yticklabels( ax(2), yticks(ax(2)) - 90); ylabel( ax(2), 'Angle [°]')

saveFigure( f, fullfile( "Z:\Emilio\SuperiorColliculusExperiments\" + ...
    "Roller\Batch7_ephys\MC\GADi52\220808_C+F_2100\ephys_E1\Figures", ...
    "Regression method" ), true, false )