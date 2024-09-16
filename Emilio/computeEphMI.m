fnOpts = {'UniformOutput', false};
tocol = @(x) x(:);
my_cat = @(x,d) cat( d, x{:} );
getMI = @(x,d) diff( x, 1, d ) ./ sum( x, d );

mi_lp = squeeze( getMI( lp_mu, 2 ) )';
popMImean_lp = [mean( mi_lp(:, tx < 0.05), 2, "omitmissing"), ...
    mean( mi_lp(:, tx >= 0.05), 2, "omitmissing")];

mi_lpu = cell( size( lPSTH ) );
uMiMatMean = zeros( sum( cellfun(@(x) size( x, 1), lPSTH ) ), 3 );
uMiMatMed = uMiMatMean;
ridx = cumsum( cellfun(@(x) size( x, 1), lPSTH ) );
r = 1;
for cs = 1:numel( lPSTH )
    mi_lpu{cs} = getMI( lPSTH{cs}, 3 );
    uMiMatMean(r:ridx(cs),:) = [repmat( cs, size( mi_lpu{cs}, 1 ), 1 ), ...
        mean( mi_lpu{cs}(:, tx < 0.05), 2, "omitmissing" ), ...
        mean( mi_lpu{cs}(:, tx >= 0.05), 2, "omitmissing" )];
    % uMiMatMed(r:ridx(cs),:) = [repmat( cs, size( mi_lpu{cs}, 1 ), 1 ), ...
    %     median( mi_lpu{cs}(:, tx < 0.05), 2, "omitmissing" ), ...
    %     median( mi_lpu{cs}(:, tx >= 0.05, "omitmissing" ), 2)];
    r = 1 + ridx(cs);
end

popMImedian = my_cat( arrayfun( @(s) ...
    median( uMiMatMean( uMiMatMean(:,1) == s, [2,3] ), 1, "omitmissing" ), ...
    unique( uMiMatMean(:,1) ), fnOpts{:} ), 1 );
popMImean = my_cat( arrayfun( @(s) ...
    mean( uMiMatMean( uMiMatMean(:,1) == s, [2,3] ), 1, "omitmissing" ), ...
    unique( uMiMatMean(:,1) ), fnOpts{:} ), 1 );

figure; 
subplot(1,5,[1,4])
boxchart( tocol( repmat( uMiMatMean(:,1), 2, 1) ), ...
    tocol( uMiMatMean(:,[2,3]) ), "Notch", "on", ...
    "GroupByColor", tocol( ones( size( uMiMatMean, 1), 1) * (1:2) ) ); 
xlim( [0, numel(lPSTH)] + 0.5 )
yline( 0, 'k--'); set( gca, "Box", "off", "Color", "none" )
subplot(1,5,5)
boxchart(popMImean, 'Notch', 'on' )
yline( 0, 'k--'); set( gca, "Box", "off", "Color", "none" )

figure; boxchart( popMImean_lp, 'Notch', 'on' );
yline( 0, 'k--'); set( gca, "Box", "off", "Color", "none" )

figure; smth = 0.015;
daviolinplot( uMiMatMean(:,[2,3]), 'groups', tocol( ones( size( uMiMatMean, 1), 1 ) * (1:2) ), ...
    'violinalpha', 0.5, 'smoothing', smth, ... 'jitter', 2, ...
    'color', [0,0.51,1] );