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

exclude_names = {'GADi13', 'GADi15', 'GADi53'};

iRN_mice = dir( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch*\MC\GADi*" );
animalFolders = arrayfun(@(f) string( expandName( f ) ), iRN_mice(:));
exclude_flags = contains( animalFolders, exclude_names );

%% Looping animals
oldMouse = "";
mc = 0; mice = [];
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
        sessOrgDirs = arrayfun(@(d) string(d.name), childFolders);
        sessOrgDirs(contains(sessOrgDirs, {'behaviour', 'ephys', ...
            'figures', 'opto'}, ctOpts{:})) = [];
        miFigDir = arrayfun(@(d) recursiveFolderSearch(expandName(d), ...
            ephFF, 'SearchType', 'expression'), childFolders, fnOpts{:}); 
        miFigDir = cat(1, miFigDir{:});
        miIdxFiles = arrayfun(@(d) dir( fullfile( d, ...
            "LogPSTH_Structure*.mat" ) ), miFigDir, fnOpts{:} ); 
        miIdxFiles = cat( 1, miIdxFiles{:} );
        if isempty(miIdxFiles)
            fprintf(1, 'No ephys analysis done! Skipping %s!\n', curDir)
            continue
        end
        % miIdxStr = arrayfun(@(x) load( expandName(x), 'logPSTH' ), miIdxFiles);
        load( expandName(miIdxFiles), 'logPSTH' );
        lp_mu = squeeze( mean( ...
            logPSTH.LogPSTH(:,:,logPSTH.indexMIComparison) ) );
        muMI = getMI( lp_mu, 2 ); muMI(isnan( muMI )) = 0;
        bmot_MI = mean( muMI( ~(logPSTH.TimeAxis < 5e-2) ) );
        sens_MI = mean( muMI( logPSTH.TimeAxis < 5e-2) );
        miVal = [sens_MI, bmot_MI];
        condNames = logPSTH.ConditionNames( logPSTH.indexMIComparison );
        condNames = join([condNames(1), "v", condNames(2)]);

        sessType = 'single';
        dataTable = table( condNames, miVal, tblOpts{:} );
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
    end
end
mice( arrayfun(@(x) isempty(x.Sessions), mice) ) = [];

behFP = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller", "MCiRNs_EphMI_sm.mat" );
svOpts = {'-mat'};
if exist(behFP, "file")
    svOpts = {'-append'};
end
%save(behFP, "mice", svOpts{:})
habFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "multi", ...
    m.Sessions), mice, fnOpts{:});
cat( 1, habFlag{:} )