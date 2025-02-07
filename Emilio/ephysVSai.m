fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);
animalPattern = '[a-zA-Z]+\d{1,}';
rsOpts = {animalPattern, 'SearchType', 'expression'};
ctOpts = {'IgnoreCase', true};
lsOpts = {'L\d+.\d+', 'match'};
ephFF = 'Ephys VW(-?\d+\.\d+)-(\d+\.\d+) RW20.00-200.00 SW(-?\d+\.\d+)-(-?\d+\.\d+)';
cond_exp = {'^Delay\s\d[.]\d{3}\s\w\s[+]\sL[0-9.]', ... iRNs
    '^Delay\s\d[.]\d{3}\s\w'}; % Continuous
cond_sel = 2;
tblOpts = {'VariableNames', {'Conditions', 'MI'}};
my_zscore = @(x, m, s) ( x - m ) ./ ( s .* (s~=0) + 1 .* (s==0) );
owfFlag = false;
m = 1e-3; k = 1e3;

cellcat = @(x,d) cat( d, x{:} );
tocol = @(x) x(:);
getMI = @(x,d) diff(x, 1, d) ./ ...
    sum(x, d).*(sum(x,d)>0) + (1.*(sum(x,d)==0 | sum(x,d)< 1e-12));
% total_var_dist = @(dmat) integral( @(x) abs( pdf( dmat(1), x ) - pdf( dmat(2), x ) ), -5, 5 );

% exclude_names = {'GADi13', 'GADi15', 'GADi53'};
exclude_names = { };
vWin = [-300, 400]*m;
% iRN_mice = dir( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch*\MC\GADi*" );
iRN_mice = dir( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch*\eOPN3\*" );
iRN_mice = iRN_mice( ~cellfun('isempty', regexp({iRN_mice.name}, '^\w{2}\d{2}' ) ) );
animalFolders = arrayfun(@(f) string( expandName( f ) ), iRN_mice(:));
exclude_flags = contains( animalFolders, exclude_names );
animalFolders = animalFolders( ~exclude_flags );
%%
% asPaths = arrayfun(@(a) dir( fullfile( a, '*F*', 'ephys*', 'cluster_info.tsv' ) ), ...
%     animalFolders, fnOpts{:} );
asPaths = arrayfun(@(a) dir( fullfile( a, '*', 'ephys*', 'cluster_info.tsv' ) ), ...
    animalFolders, fnOpts{:} );
bad_flag = cellfun(@(x) contains( {x.folder}, 'bad' ), asPaths, fnOpts{:} );
notF_flag = cellfun(@(x) contains( {x.folder}, 'f' ), asPaths, fnOpts{:} );
asPaths = cellfun(@(a,f,w) a(~f&~w), asPaths, bad_flag, notF_flag, fnOpts{:} );
empty_flag = cellfun('isempty', asPaths );
asPaths = asPaths(~empty_flag); animalFolders = animalFolders(~empty_flag);
clInfo = cellfun(@(a) arrayfun(@(s) getClusterInfo( expandName( s ) ), ...
    a, fnOpts{:} ), asPaths, fnOpts{:} );

%%
Nu = cellfun(@(x) cellfun(@(y) sum( y.ActiveUnit ), x ), clInfo, fnOpts{:} );
low_unit_flag = cellfun(@(x) x < 10, Nu, fnOpts{:} );
low_unit_animal = cellfun(@(x) all(x < 10), Nu );
asPaths = cellfun(@(a,f) a(~f), asPaths, low_unit_flag , fnOpts{:} );
animalFolders = animalFolders( ~low_unit_animal );
clInfo =  cellfun(@(a,f) a(~f), clInfo, low_unit_flag , fnOpts{:} );
Nm = numel( clInfo );
Nu =  cellfun(@(a,f) a(~f), Nu, low_unit_flag , fnOpts{:} );
% clInfo = arrayfun(@(x) getClusterInfo( expandName( ...
%     dir( fullfile( x, '*F*', 'ephys*', 'cluster_info.tsv' ) ) ) ), animalFolders, ...
%     fnOpts{:});
%% Pooling ephys and behaviour
% Initialising
Nspm = cellfun("size", clInfo, 1 ); % Number of sessions per mouse
Nexp = sum( Nspm );
Nu_all = cat( 1, Nu{:} );
Nuinit = cumsum( [1; Nu_all(1:end-1) ] );
Nuend = cumsum( Nu_all );
PSTHall_mu = zeros( 700, sum( Nu_all ), 2 );
brAll = [];
PSTHall = cell( Nexp, 2 );
uSig = cell( Nexp, 1 ); uID = uSig;
uMod = uSig;
uMI = uSig; ce = 1;
%%
rstVars2load = {'relativeSpkTmsStruct', 'configStructure'};
afVars2load = {'Conditions', 'fs'};
for ca = 1:Nm
    for cs = 1:Nspm(ca)
        data_dir = string( getParentDir( asPaths{ca}(cs).folder, 1 ) );
        % condStruct = load( expandName( dir( fullfile( data_dir, "*\*analysis.mat" ) ) ), afVars2load{:} );
        rstPath = dir( fullfile( data_dir, "*", "*RW20.00-200.00*(unfiltered) RelSpkTms.mat" ) );
        brPath = dir( fullfile( data_dir, "*","BehaviourResult*.mat" ) );
        mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Res VW* RW20.00-200.00 ms SW*PuffAll.mat") );
        brVars2load = 'behRes';
        if isempty( brPath )
            brPath = dir( fullfile( data_dir, "*\Simple summary.mat" ) );
            brVars2load = 'summStruct';
        end
        brStruct = load( expandName( brPath ), brVars2load );
        behRes = brStruct.(brVars2load);
        ctrlSub = ismember( string( {behRes.ConditionName} ), 'Control Puff' );
        % cf_flags = regexp( string( {behRes.ConditionName} ), {'Control Puff', ...
        %     'Delay\s\d[.]\d{3}\s\w\s[+]\sL[0-9.]'} );
        % cf_flags = contains( string( {behRes.ConditionName} ), 'Control Puff');
        cf_flags = regexp( string( {behRes.ConditionName} ), ...
             cond_exp{cond_sel} );
        cf_flags = ctrlSub | ~cellfun( 'isempty', cf_flags );

        brAll = cat( 1, brAll, behRes(cf_flags) );
        if ~isempty(rstPath) || numel( rstPath ) == 1
            rstCont = load( expandName( rstPath ), rstVars2load{:} );
        else
            fprintf(1, 'Either empty or more than 1 file found!\n');
            disp( {rstPath.name} )
            continue
        end
        mftype = 1;
        mfVars2load = {'Results', 'gclID', 'Counts', 'configStructure'};
        if isempty( mfPath )
            fprintf( 1, 'Res file not found!\n')
            mfPath = dir( fullfile( data_dir, "ephys*", "Results", "Map*.mat") );
            mfVars2load = {'keyCell', 'resMap'};
            if isempty( mfPath )
                fprintf( 1, 'Unable to load unit response information!\n')
                disp( mfPath )
                continue
            end
            mftype = 2;
        end
        mfStruct = load( expandName( mfPath ), mfVars2load{:} );
        if mftype == 1
            Results = mfStruct.Results;
            % Taking only control and laser frequency
            cnfStruct = mfStruct.configStructure;
            Nccond = numel( cnfStruct.ConsideredConditions );
            valid_cond_subs = 1:Nccond;
            ctrlSub = ismember( cnfStruct.ConsideredConditions, 'Control Puff' );
            cf_flags = regexp( cnfStruct.ConsideredConditions, ...
                cond_exp{cond_sel} );
            cf_flags = ctrlSub | ~cellfun( 'isempty', cf_flags );
            valid_cond_subs = valid_cond_subs( cf_flags );
            cmbSubs = cellcat(arrayfun(@(x) sscanf( x.Combination, '%d %d'), ...
                Results , fnOpts{:} ), 2 );
            configFlag = cmbSubs(1,:) == cmbSubs(2,:) & ...
                any(cmbSubs(1,:) == valid_cond_subs(:), 1);
            cond_ordr = cmbSubs( 1, configFlag );
            gclID = mfStruct.gclID;

            Counts = mfStruct.Counts;
            Counts = arrayfun(@(c) squeeze( mean( cellcat( Counts(c,:), 3 ), 2 ) ), ...
                valid_cond_subs, fnOpts{:} );
            % configFlag = cellfun(@(x) ~isempty(x), regexp( {Results.Combination}, ...
            %     '1\s1\ssignrank', 'ignorecase' ) );
            sig_aux = arrayfun(@(x) x.Activity(1).Pvalues, Results(configFlag), fnOpts{:} );
            mod_aux = arrayfun(@(x) x.Activity(1).Direction, Results(configFlag), fnOpts{:} );
            mi_aux = cellfun(@(x) getMI( x, 2 ), Counts, fnOpts{:} );
            uSig{ce} = cellcat( sig_aux(cond_ordr), 2);
            uMod{ce} = cellcat( mod_aux(cond_ordr), 2);
            uMI{ce} = cellcat( mod_aux(cond_ordr), 2 );
            uID{ce} = gclID;
        else
            resMap = mfStruct.resMap; keyCell = mfStruct.keyCell;
            configFlag = contains( keyCell(:,1), 'RW20.00-200.00' ) & ...
                contains( keyCell(:,4), 'Control Puff' );
            if sum( configFlag ) ~= 1
                fprintf( 1, 'Cannot process this keycell... \n')
                disp( keyCell )
            end
            uSig{ce} = resMap( keyCell{configFlag,:} );
            uMod{ce} = zeros( Nu(ca), 1 ); uMI{ca} = uMod{ca};
            uID{ce} = clInfo{ca}{clInfo{ca}.ActiveUnit==1, "cluster_id" };
        end
        rstStruct = rstCont.relativeSpkTmsStruct;
        confStruct = rstCont.configStructure;
        % Conditions = condStruct.Conditions;
        % fs = condStruct.fs;
        if any( confStruct.Viewing_window_s ~= vWin )
            confStruct.Viewing_window_s = vWin;
        end
        [PSTH, trial_tx, Na] = getPSTH_perU_perT( ...
            rstStruct(valid_cond_subs), confStruct );
        a = Nuinit(ce); b = Nuend(ce);
        idx = a:b;
        PSTHall_mu(:,idx,:) = cellcat( cellfun(@(x) squeeze( mean( x, 1 ) ), ...
            PSTH, fnOpts{:} ), 3 );
        PSTHall(ce,:) = PSTH;
        ce = ce + 1;
    end
end
clearvars PSTH brStruct ctrlSub rstStruct confStruct brPath brVars2load ...
    rstCont a b idx;
ai_pt = arrayfun(@(s) getAIperTrial( s ), brAll, fnOpts{:} );

%% Sliding window analysis for behaviour correlation
% Idea is to slide a time window per unit per trial for getting an R²
slid_win_length = 20*m; time_slide = 5*m;
time_init = -50*m; time_stop = 400*m;
Nrs = (time_stop - time_init - slid_win_length) / time_slide;
r_squared = cell( Nexp, 1 );
Nt = cellfun( "size", PSTHall, 1 );
parfor cexp = 1:Nexp
    r_squared{cexp} =  zeros( Nu_all(cexp), Nrs, 2 );
    for cu = 1:Nu_all(cexp)
        aux_rs = zeros(1, Nrs, 2 );
        for ct = 1:2
            cw = time_init + [0, slid_win_length];
            ci = 1;
            aux_trial = cellcat( arrayfun(@(t) conv( PSTHall{cexp,ct}(t, :, cu ), ...
                gausswin( 5 ), "same" ), 1:Nt(cexp,ct), 'UniformOutput', false ), 1 );
            while cw(2) <= time_stop
                act_mu = mean( aux_trial(:, my_xor( trial_tx < cw ) ) , 2 );
                if ( sum( act_mu == 0 ) / numel(act_mu ) ) < 0.4
                    aux_mdl = fitlm( zscore( act_mu )', zscore( ai_pt{cexp,ct} )', 'poly1' );
                    aux_rs(1,ci,ct) = aux_mdl.Rsquared.Ordinary;
                end
                cw = cw + time_slide; ci = ci + 1;
            end
        end
        r_squared{cexp}(cu,:,:) = aux_rs;
    end
end
% r_squared_cat = cat( 1, r_squared{:} );
%% Time resolved boxplots for all experiments
td = 30;
slid_win_length = td*m; time_slide = td*m;
time_init = -160*m; time_stop = 400*m;
Nrs = floor( (time_stop - time_init - slid_win_length) / time_slide );
r_squared_pexp = zeros( Nexp, Nrs, 2 );
parfor cexp = 1:Nexp
    for ct = 1:2
        cw = time_init + [0, slid_win_length];
        for ci = 1:Nrs
            act_mu = mean( PSTHall{cexp,ct}(:, my_xor( trial_tx < cw ),: ) , [2,3] );
            aux_mdl = fitlm( zscore( act_mu )', zscore( ai_pt{cexp,ct} )', 'poly1' );
            r_squared_pexp(cexp,ci,ct) = aux_mdl.Rsquared.Ordinary;
            cw = cw + time_slide;
        end
    end
end
aux_mdl = fit_poly( [1,Nrs], [time_init, time_init + Nrs*slid_win_length] ...
    + (slid_win_length/2)*[1,-1], 1 );
b_tx = ( ( 1:Nrs )'.^[1,0] ) * aux_mdl;
% r_squared_cat = cat( 1, r_squared_pexp{:} );
%% 
bxOpts = {'JitterOutliers', 'on', 'MarkerStyle', '.', 'MarkerColor', 'k',...
    'BoxFaceColor', 'k', 'BoxWidth', k*(time_slide)/2, ...
    'Notch', 'off' };
ttl = ["Laser OFF", "Laser ON"];
f = figure('Color', 'w'); t = createtiles(f, 2, 2); 
ax = gobjects( 3, 1 );
for ct = 1:2
    ax(ct) = nexttile( t );
    boxchart(ax(ct), tocol( ones(Nexp,1)*b_tx' * k ), ...
        tocol(r_squared_pexp(:,:,ct)), bxOpts{:} )
    hold( ax(ct), 'on');
    line(ax(ct), k*b_tx, median( r_squared_pexp(:,:,ct), 1 ), 'Color', 'k', ...
        'LineWidth', 2 )
    xline( ax(ct), [0,50,200], 'r--')
    xlabel( ax(ct), 'Time [ms]' )
    xlim( ax(ct), k*(b_tx([1,end]) + [-1;1]*slid_win_length/2) )
    if ct ~=2
        ylabel( ax(ct), 'R² per window' )
    else
        set( get( ax(ct), 'YAxis'), 'Visible', 'off' )
    end
    cleanAxis( ax(ct) );
    title( ax(ct), ttl(ct) )
    
end
ct = ct + 1;
title( t, sprintf('Time-resolved_{%d ms} R² population (per experiment, eOPN3)', td ) )
ax(ct) = nexttile( t, [1,2] );
lObjs = line( ax(ct), b_tx*k, squeeze( median( r_squared_pexp, 1 ) ) );
legend( lObjs, ttl, 'Color', 'none', 'Box', 'off', 'Location', 'best');
xlabel( ax(ct), 'Time [ms]' )
ylabel( ax(ct), 'R²' )
ytickangle( ax, 90 )
set( ax, 'TickDir', 'out' )
linkaxes(ax, 'x'); linkaxes( ax(1:2), 'y' )

%%
% Comparing conditions agains each other 
p_cond_per_window = arrayfun(@(x) signrank( ...
    squeeze( r_squared_pexp(:,x,1) ), ...
    squeeze( r_squared_pexp(:,x,2) ) ), ...
    1:Nrs );
p_th = [0.05; 0.01; 0.001];
astk = sum( p_cond_per_window < p_th );
text( ax(ct), b_tx*k, max( median( r_squared_pexp ), [], 3 ) * 1.15, ...
    cellfun(@(x) join(x), arrayfun(@(x) repmat( "\ast", 1, x ), astk, fnOpts{:} ) ), ...
    "HorizontalAlignment", "center", "VerticalAlignment", "middle", "FontSize", 17 )
%%  Statistics on the time-resolved R²
tbl2 = cell(2, 1 );
for ct = 1:2
    [p, tbl, stats] = kruskalwallis( squeeze(r_squared_pexp(:,:,ct) ), ...
        string( b_tx(:) * k ) );
    figure; tbl2{ct} = multcompare( stats );
    sum( tbl2{ct}(:,end) < 0.05 )
end
