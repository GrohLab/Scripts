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
my_zscore = @(x, m, s) ( x - m ) ./ ( s .* (s~=0) + 1 .* (s==0) );
getMI = @(x,d) diff(x, 1, d)./sum(x, d);
total_var_dist = @(dmat) integral( @(x) abs( pdf( dmat(1), x ) - pdf( dmat(2), x ) ), -5, 5 );
tocol = @(x) x(:);
m = 1e-3;
exclude_names = {'GADi13', 'GADi15', 'GADi53'};

params = struct( 'relative_window', [-1,1]*800*m, 'delay_window', ...
    [-1,1]*100*m, 'bin_size', 5*m, 'kfold', 20 );

pc = parcluster('local');
if ~strcmp( computer, 'PCWIN64')
    home_path = '/gpfs/bwfor/home/hd/hd_hd/hd_bf154/';
    repo_paths = cellfun(@(x) char( fullfile( home_path, x) ), ...
        {'NeuroNetzAnalysis', 'AuxiliaryFuncs', 'Scripts'},  fnOpts{:} );
    addpath( repo_paths{:} )
    roller_path = "/mnt/sds-hd/sd19b001/Emilio/SuperiorColliculusExperiments/Roller";
else
    roller_path = "Z:\Emilio\SuperiorColliculusExperiments\Roller";
end
parpool( pc )

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

        [mdl, params, DX] = regressEphysVSBehaviour( data_path, params );
        if isempty(DX)
            continue
        end
        mdl_mu = squeeze( mean( mdl, 2 ) );
        y_trials = reshape( DX{1}, params.Nb, params.Nr, params.Ns );
        y_pred = DX{2} * mdl_mu;
        y_ptrials = reshape( y_pred, params.Nb, params.Nr, params.Ns );

        SSEt = squeeze( sum( ( y_trials - y_ptrials ).^2, 1 ) );
        SSTt = squeeze( sum( ( y_trials - mean( y_trials, 1 ) ).^2, 1 ) );
        r_sq_trials = 1 - (SSEt./SSTt);

        SSE = sum( (DX{1} - y_pred).^2 );
        SST = sum( (DX{1} - mean( DX{1}, 1 ) ).^2 );
        r_sq = 1 - ( SSE./ SST );

        dataTable = table( r_sq, {r_sq_trials}, ...
            'VariableNames', {'R_squared', 'R_squared_trials'});
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

behFP = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller", "MCiRNs_regression_sm.mat" );
svOpts = {'-mat'};
if exist(behFP, "file")
    svOpts = {'-append'};
end
% save(behFP, "mice", svOpts{:})
habFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "multi", ...
    m.Sessions), mice, fnOpts{:});
cat( 1, habFlag{:} )
