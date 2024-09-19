
%{
res_gof = zeros( 15, 1 );
% feps = zeros( 15, 1 );
parfor ii = 1:15
    trainIdx = training( cvpart, ii );
    testIdx = test( cvpart, ii );
    Xtrain = gpuArray( X( any( tr_ID == diag( trainIdx * (1:Nr) )', 2 ), : ) );
    Ytrain = gpuArray( y( any( tr_ID == diag( trainIdx * (1:Nr) )', 2 ), 1 ) );
    Xtest = gpuArray( X( any( tr_ID == diag( testIdx * (1:Nr) )', 2 ), : ) );
    Ytest = gpuArray( y( any( tr_ID == diag( testIdx * (1:Nr) )', 2 ), 1 ) );

    mdl = fitglm( Xtrain, Ytrain, 'linear', 'Distribution', 'Normal' );
    y_pred = feval( mdl, Xtest );
    res_gof(ii) = goodnessOfFit( Ytest, y_pred, 'MSE' );
    % feps(ii) = loss( mdl, Xtest, Ytest );
end
%}
%%
fnOpts = {'UniformOutput', false};
tocol = @(x) x(:);
roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
data_path = fullfile( roller_path, "Batch18_ephys\MC\GADi43\240227_C+F_2200" );
eph_path = fullfile( data_path, "ephys_E1" );
beh_path = fullfile( data_path, "Behaviour" );
load( fullfile( beh_path, "BehaviourSignals2024-02-27T11_26_07+T11_43_01.mat" ) )
load( fullfile( beh_path, "RollerSpeed2024-02-27T11_26_07+T11_43_01.mat" ) )
% load( fullfile( eph_path, "GADi43_C+F_2200_all_channels.mat" ) )
load( fullfile( eph_path, ...
    "GADi43_C+F_2200 RW20.00-50.00 SW-180.00--150.00 VW-300.00-400.00 ms PuffAll (unfiltered) RelSpkTms.mat" ), "configStructure" )
load( fullfile( eph_path, "GADi43_C+F_2200analysis.mat" ) )
load( fullfile( eph_path, "GADi43_C+F_2200_Spike_Times.mat" ) )
stop_time = length( Triggers.Whisker )/ fs;
%%
behSignals = [behDLCSignals, vf];
mdl_btx = fit_poly( [1, size( behSignals, 1 )], [0, size( behSignals, 1 )/fr] + [1,-1] * (1/fr), 1 );
btx = (1:size( behSignals, 1 ))'.^[1,0] * mdl_btx;
my_xor = @(x) xor( x(:,1), x(:,2) );
my_cat = @(x,d) cat( d, x{:} );

m = 1e-3;
% time_limits = [0, length(behSignals)/fr];
cS = configStructure;
rel_win = [-1, 1]*0.8;
del_win = [-50, 50]*m;
bin_size = 5*m;
cS.BinSize_s = bin_size;
Nb = ceil( diff( rel_win )/ bin_size );
% Nb = ceil( diff( time_limits ) / bin_size );
Nu = numel( spike_times );
% cons_time = my_xor( btx > time_limits );
Ns = size( behSignals, 2 );
%%
bin_edges = 0:bin_size:stop_time;
bin_centres = mean( [bin_edges(1:end-1); bin_edges(2:end)] );
Ntb = length( bin_centres );
hstOpts = {'Normalization', 'countdensity'};
binned_spikes = cellfun(@(s) histcounts( s, bin_edges, hstOpts{:}), ...
    spike_times, fnOpts{:} );
binned_spikes = cat( 1, binned_spikes{:} );

binned_beh = zeros( Ntb, Ns );
parfor b = 1:Ntb
    idx = my_xor( btx(:) < bin_edges(b:b+1) );
    binned_beh(b,:) = mean( behSignals( idx , : ), 1 );
end

%% Design matrix for a set of trials
time_limits = Conditions(3).Triggers(:,1)./fs + rel_win;
Nr = size( time_limits, 1 );
Nd = ceil( diff( del_win ) / bin_size );
auX = zeros( Nb*Nr, Nu, Nd );

cwin = arrayfun(@(x) linspace( time_limits(x,1) + (bin_size/2), ...
    time_limits(x,2) - (bin_size/2), Nb )', (1:Nr)', fnOpts{:} );
cwin = cat( 1, cwin{:} );

bin_ax = cwin + linspace( del_win(1)+(bin_size/2), ...
    del_win(2)-(bin_size/2), Nd );
% tr_ID = ceil( ( 1:(Nr*Nb) )' / Nb );
parfor r = 1:(Nr*Nb)
    tempC = my_cat( arrayfun( @(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax(r,:) ), 1:Nu, fnOpts{:} ), 1);
    auX( r, :, :) = tempC;
end

X = reshape( auX, [], Nu*Nd ); clearvars auX;
X2 = [ ones( Nb*Nr, 1), X];
Xp = X2;
%% Multivariate regression response matrix
%X2 = [ ones( Nb*Nr, 1), X];
%lmObjs = cell( Ns, 1 );
y = zeros( Nb*Nr, Ns);
for r = 1:Nr
    idx = (r-1)*Nb + (1:Nb);
    aux = arrayfun(@(s) interp1( bin_centres, binned_beh(:,s), ...
        (1:Nb)'*bin_size + time_limits(r,1) ), (1:Ns), fnOpts{:} );
    aux = cat( 2, aux{:} );
    y(idx,:) = aux;
end
%%
tr_ID = tocol( ones( Nb, 1 ) * (1:Nr) );
rmse = zeros( 15 , 1 ); mdl = cell( size( rmse ) );
Nk = round( Nr*0.15 ); 
parfor ii = 1:Nk
    testTrials = sort( randperm( Nr, Nk ) );
    trainingTrials = setdiff( 1:Nr, testTrials );
    trainingIdx = any( tr_ID == trainingTrials(:)', 2 );
    testIdx = ~trainingIdx;

    mdl{ii} = fitlm( X(trainingIdx,:), y(trainingIdx,1) );
    y_pred = predict( mdl{ii}, X(testIdx,:) );
    rmse(ii) = sqrt( mean( ( y(testIdx,1) - y_pred ).^2 ) );
end

[~, min_error] = min(rmse);
y_all_pred = predict( mdl{min_error}, X );

y_1 = reshape( y(:,1), Nb, Nr );
y_all_pred = reshape( y_all_pred, Nb, Nr );

%% 
cvk = 15;
Nk = round( Nr*0.15 );
rmse = zeros( cvk , Ns ); 
mdl = zeros( size( X, 2 )+Ns, size( y, 2 ), cvk );
idxs = zeros( cvk, Nk ); %[zy, y_mu, y_sig] = zscore(y, 0, 1);
X2 = [[eye(Ns); zeros( size(X,1) - Ns, Ns )], X];
for ii = 1:cvk
    fprintf(1, 'K:%d\n', ii)
    testTrials = sort( randperm( Nr, Nk ) );
    idxs(ii,:) = testTrials;
    trainingTrials = setdiff( 1:Nr, testTrials );
    trainingIdx = any( tr_ID == trainingTrials, 2 );
    testIdx = ~trainingIdx;

    mdl(:,:,ii) = mvregress( gpuArray( X2(trainingIdx,:) ), ...
        gpuArray( y(trainingIdx,:) ) );
    y_pred = X2(testIdx,:) * mdl(:,:,ii);
    rmse(ii,:) = sqrt( mean( ( y(testIdx,:) - y_pred ).^2 ) );
end

[~, min_error] = min(rmse,[],1);
y_all_pred = X2 * squeeze( mean( mdl, 3 ) );

y_trials = reshape( y, Nb, Nr, Ns );
y_all_pred = reshape( y_all_pred, Nb, Nr, Ns );
%% Training

ho_trials = randperm( Nr, round( Nr*0.1 ) );
testIdx = any( tr_ID == sort(ho_trials), 2 );
cv_kf = cvpartition( tr_ID( ~testIdx ), "KFold", 15 );

Xtrain = X(~testIdx, :); ytrain = y(~testIdx,1);

[mdl, fitInfo] = lassoglm( X(~testIdx,:), y(~testIdx,1), 'normal', ...
    'CV', cv_kf, 'Lambda', logspace( -5, 3, 64 ), ...
    'Options', statset('UseParallel', true ), ...
    'Alpha', eps );

w_vec = [fitInfo.Intercept(fitInfo.IndexMinDeviance);
    mdl(:,fitInfo.IndexMinDeviance)];
y_pred = glmval( w_vec, X, "identity" );
y_pred = reshape( y_pred, Nb, Nr );
clrMap = [0.15*ones(1,3); 0.85,0.51,0.15 ];
figure; lObj = line( 1:(Nb*numel( ho_trials )), [y(testIdx,1), y_pred(:)] );
arrayfun(@(ii,x) set( x, 'Color', clrMap(ii,:) ), (1:numel(lObj))', lObj(:) )

createtiles = @(f,nr,nc) tiledlayout( f, nr, nc, ...
    'TileSpacing', 'Compact', 'Padding', 'tight');
cleanAxis = @(x) set( x, "Box", "off", "Color", "none" );

fig = figure( 'Color', 'w' ); t = createtiles( fig, 10, 2);
nexttile([8, 1]);
imagesc( ttx, [], y_1' - 90); colormap(inferno)
xline(0, 'LineStyle', '--', 'Color', 0.85*ones(1,3) )
nexttile([8, 1]);
imagesc( ttx, [], y_pred' - 90); colormap(inferno)
set( get( gca, 'YAxis' ), 'Visible', 'off' )
nexttile(t);
line( ttx, mean( y_1, 2 ) - 90, 'Color', 0.15*ones(1,3), 'LineWidth', 1.5 )
nexttile(t);
line( ttx, mean( y_pred, 2 ) - 90, 'Color', [0.85, 0.51, 0.15] , 'LineWidth', 1.5 )
axs = get( t, "Children" );
linkaxes( axs, 'x')
xlim(ttx([1,end]))
arrayfun(@(x) set( get( x, "XAxis" ), "Visible", "off" ), axs(3:4) )
arrayfun(@(x) xticklabels( x, xticks(x) / m ), axs(1:2) )

rmse = mean( ( y_1 - y_pred ).^2, 1 );


wtx = (del_win(1) + bin_size/2):bin_size:(del_win(2) - bin_size/2);
ttx = (rel_win(1) + bin_size/2):bin_size:(rel_win(2) - bin_size/2);

figure; imagesc( wtx, [], reshape( w_vec(2:end), Nd, Nu )' )

%% Design matrix for the whole experiment
Nb = ceil( diff( bin_edges([1,end]) )/ bin_size );
auX = zeros( Nb, Nu, Nd );
parfor b = 1:Nb
    cwin = bin_centres(b) + del_win;
    bin_ax = linspace( cwin(1), cwin(2), Nd );
    tempC = arrayfun(@(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax ), 1:Nu, fnOpts{:} );
    tempC = cat( 1, tempC{:} );
    tempC( isnan(tempC) ) = 0;
    auX( b, :, :) = tempC;
end
X = reshape( auX, [], Nu*Nd );
X2 = [ ones( Nb*Nr, 1), X];
Xa = X2;
