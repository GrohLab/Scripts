%#ok<*AGROW,*SAGROW>
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);
animalPattern = '[a-zA-Z]+\d{1,}';
rsOpts = {animalPattern, 'SearchType', 'expression'};
ctOpts = {'IgnoreCase', true};
lsOpts = {'L\d+.\d+', 'match'};
ephFF = 'Ephys VW(-?\d+\.\d+)-(\d+\.\d+) RW20.00-200.00 SW(-?\d+\.\d+)-(-?\d+\.\d+)';
tblOpts = {'VariableNames', {'Conditions', 'MI'}};
% my_zscore = @(x, m, s) ( x - m ) ./ ( s .* (s~=0) + 1 .* (s==0) );
% getMI = @(x,d) diff(x, 1, d)./sum(x, d);
% total_var_dist = @(dmat) integral( @(x) abs( pdf( dmat(1), x ) - pdf( dmat(2), x ) ), -5, 5 );
getRMSE = @( r, x, d ) sqrt( mean( ( r - x ).^2, d, "omitmissing" ) );
tocol = @(x) x(:);
m = 1e-3;
exclude_names = {'GADi13', 'GADi15', 'GADi53'};

params = struct( 'relative_window', [-1,1]*800*m, 'delay_window', ...
    [-1,1]*100*m, 'bin_size', 5*m, 'kfold', 20 );
Nbins = diff(params.relative_window)/params.bin_size;
time_mdl = fit_poly( [1,Nbins], params.relative_window + ...
    [1,-1]*(params.bin_size/2), 1 );
tx = ( ( 1:Nbins)'.^[1,0] ) * time_mdl;
sponFlags = tx < 0;
getTimeBoAT = @(s,f,p) reshape( s(f,:,:), size(s,2) * sum( f ), p.Ns );

pc = parcluster('local');
if ~strcmp( computer, 'PCWIN64')
    home_path = '/gpfs/bwfor/home/hd/hd_hd/hd_bf154/';
    repo_paths = cellfun(@(x) char( fullfile( home_path, x) ), ...
        {'NeuroNetzAnalysis', 'AuxiliaryFuncs', 'Scripts'},  fnOpts{:} );
    cellfun(@(x) addpath( genpath( x ) ), repo_paths )
    roller_path = "/mnt/sds-hd/sd19b001/Emilio/SuperiorColliculusExperiments/Roller";
else
    roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
end
try
    parpool( pc );
catch
end

iRN_mice = dir( fullfile( roller_path, "Batch*", "MC", "GADi*" ) );
animalFolders = arrayfun(@(f) string( expandName( f ) ), iRN_mice(:));
exclude_flags = contains( animalFolders, exclude_names );
%% Looping animals
oldMouse = "";
mc = 0; mice = []; lp_mu = []; lPSTH = [];
sessType = 'single';
for cad = tocol(animalFolders(~exclude_flags))'
    [structPath, currMouse] = fileparts(cad);
    [~, structName] = fileparts(structPath);
    if string(oldMouse) ~= string(currMouse)
        oldMouse = currMouse;
        mice = [mice; struct('Name', currMouse, 'Sessions',[], ...
            'Structure', structName)];
        mc = mc + 1;
        sc = 0; oldSess = ""; oldDepth = "";
    end
    fprintf(1, 'Mouse %s ', currMouse)
    sessDirs = getSubFolds(cad);
    % Just date sessions
    onlyDateSessFlag = arrayfun(@(x) string(regexp(x.name, '[0-9]{6}', ...
        'match')), sessDirs, fnOpts{:});
    sessDirs(cellfun(@isempty, onlyDateSessFlag)) = [];
    for csd = sessDirs(:)'
        curDir = expandName(csd);
        sessDateDepth = regexp(csd.name, '(\d{6}).*_(\d{4})?', 'tokens', 'once');
        if ~isempty( sessDateDepth )
            currSess = sessDateDepth{1};
            if isempty( sessDateDepth{2} )
                depthSess = '';
            else
                depthSess = sessDateDepth{2};
            end
        else
            currSess = regexp(csd.name, '(\d{6})', 'tokens', 'once');
            depthSess = '';
            if isempty(currSess)
                fprintf( 1, "Unable to get session date and depth\n" );
                fprintf( 1, "Skipping: %s %s\n", currMouse, csd.name )
                continue
            end
        end
        childFolders = getSubFolds(curDir);
        sessOrgDirs = arrayfun(@(d) string(d.name), childFolders );
        sessOrgDirs( ~contains(sessOrgDirs, {'behaviour', 'ephys', ...
            'figures', 'opto'}, ctOpts{:}) ) = [];
        if isempty(sessOrgDirs)
            continue
        end
        data_path = curDir;
        fprintf(1, ', Session %s\n', currSess )
        try
            [mdl, params, DX] = regressEphysVSBehaviour( data_path, params );
        catch ME
            display(ME.message)
            continue
        end
        if isempty(DX) || (sum( isnan(mdl), "all" ) / numel(mdl) ) > 0.05 || ...
                any( cellfun(@isempty, DX) )
            continue
        end
        %%
        mdl_mu = squeeze( mean( mdl, 2 ) );
        y_trials = reshape( DX{1}, params.Nb, params.Nr, params.Ns );
        y_pred = DX{2} * mdl_mu;
        y_ptrials = reshape( y_pred, params.Nb, params.Nr, params.Ns );
        
        r_sq_trials = goodnessFit2( y_trials, y_ptrials, 1 );
        r_sq = goodnessFit2( DX{1}, y_pred, 1 );

        y_trials_pre = getTimeBoAT(y_trials, sponFlags, params);
        y_ptrials_pre = getTimeBoAT(y_ptrials, sponFlags, params);
        r_sq_pre = goodnessFit2(y_trials_pre, y_ptrials_pre, 1);

        y_trials_post = getTimeBoAT( y_trials, ~sponFlags, params );
        y_ptrials_post = getTimeBoAT( y_ptrials, ~sponFlags, params );
        r_sq_post = goodnessFit2(y_trials_post, y_ptrials_post, 1);

        y_lpred = DX{3} * mdl_mu;
        y_lptrials = reshape( y_lpred, params.Nb, [], params.Ns );
        y_ltrials = reshape( DX{4}, size( y_lptrials ) );
        
        r_sq_l = goodnessFit2( DX{4}, y_lpred, 1 );

        y_trials_pre = getTimeBoAT(y_ltrials, sponFlags, params);
        y_ptrials_pre = getTimeBoAT(y_lptrials, sponFlags, params);
        r_sq_lpre = goodnessFit2(y_trials_pre, y_ptrials_pre, 1);

        y_trials_post = getTimeBoAT( y_ltrials, ~sponFlags, params );
        y_ptrials_post = getTimeBoAT( y_lptrials, ~sponFlags, params );
        r_sq_lpost = goodnessFit2(y_trials_post, y_ptrials_post, 1);
        
        rmse_laser = getRMSE( DX{4}, y_lpred, 1 );

        dataTable = table( {[r_sq;r_sq_pre;r_sq_post],[r_sq_l;r_sq_lpre;r_sq_post]}, ...
            {r_sq_trials}, {params.fit_error}, ...
            rmse_laser, 'VariableNames', {'R_2_p_L', 'R_squared_trials', ...
            'RMSE_c', 'RMSE_l'} );
        if ( string(oldSess) ~= string(currSess) ) || ...
                ( string(oldDepth) ~= string(depthSess) )
            oldSess = currSess;
            oldDepth = depthSess;
            auxStruct = struct('Date', currSess, ...
                'DataTable', dataTable, 'Type', sessType, ...
                'Depth', depthSess);
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
        close all
    end
end
mice( arrayfun(@(x) isempty(x.Sessions), mice) ) = [];

behFP = fullfile( roller_path, "MCiRNs_reconstruction_sm.mat" );
svOpts = {'-mat'};
if exist(behFP, "file")
    svOpts = {'-append'};
end
save(behFP, "mice", svOpts{:})
habFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "multi", ...
    m.Sessions), mice, fnOpts{:});
cat( 1, habFlag{:} )
