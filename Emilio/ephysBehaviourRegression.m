
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

% time_limits = [0, length(behSignals)/fr];
cS = configStructure;
rel_win = [-0.8, 0.8];
del_win = [-50, 50]*1e-3;
bin_size = 5e-3;
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
    tempC = arrayfun(@(u) interp1( bin_centres, binned_spikes(u,:), ...
        bin_ax(r,:) ), 1:Nu, fnOpts{:} );
    tempC = cat( 1, tempC{:} );
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

tr_ID = tocol( ones( Nb, 1) * (1:Nr) );
cv_kf = cvpartition( tr_ID, "KFold", 15 );

[mdl, fitInfo] = lassoglm( X, Y, 'normal', ...
    'CV', cv_kf, 'Lambda', logspace( -5, 3, 64 ), ...
    'Options', statset('UseParallel', true ), ...
    'Alpha', eps );

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
