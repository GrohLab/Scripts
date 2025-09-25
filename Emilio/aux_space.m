aux_mat = zeros( sum( any( ~isnan( PTX_mat(:,[2,3],:) ), 1 ), 'all' ), 2 );
cm = 1;
for cm2 = 1:size( PTX_mat, 3 )
    dat_flag = any( ~isnan( PTX_mat(:,[2,3],cm2) ), 1 );
    if any( dat_flag )
        dat_flag2 = any( ~isnan( PTX_mat(:,[2,3],cm2) ), 2 ) ;
        aux_mat(cm,:) = PTX_mat(dat_flag2,[true, dat_flag],cm2);
        cm=cm+1;
    end
end
%%
fowFlag = false;
batch_dir = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch6_beh+Muscimol" );
exp_dirs = dir( fullfile( batch_dir, "Musc", "*", "*" ) );
exp_dirs( any( string( {exp_dirs.name} ) == ["."; ".."], 1 ) ) = [];
get_full_path = @(x) fullfile( x.folder, x.name);
m = 1e-3;
for cep = 1:numel(exp_dirs)
    exp_path = get_full_path( exp_dirs(cep) );
    eph_path = dir( fullfile( exp_path, "ephys*" ) );
    eph_path( ~[eph_path.isdir] ) = [];
    beh_path = fullfile( exp_path, "Behaviour" );

    if ~isempty( eph_path )
        eph_path = get_full_path( eph_path );
        figure_path = fullfile( eph_path, "Figures" );
        af_path = dir( fullfile( eph_path , "*analysis.mat" ) );
        af_path = get_full_path( af_path );
    elseif exist( beh_path, "dir" )
        af_path = get_full_path( dir( fullfile( beh_path, "*analysis.mat") ) );
        figure_path = fullfile( beh_path, "Figures" );
    else
        beh_path = exp_path;
        af_path = get_full_path( dir( fullfile( beh_path, "*analysis.mat") ) );
        figure_path = fullfile( beh_path, "Figures" );
    end

    [~, af_name] = fileparts( af_path );
    expName = extractBefore(af_name, "analysis");
    load( af_path, "Conditions", "fs")

    fnOpts = {'UniformOutput', false};
    axOpts = {'Box','off','Color','none'};
    lgOpts = cat( 2, axOpts{1:2}, {'Location','best'} );

    ldFlag = false;
    try
        load( get_full_path( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
    catch
        ldFlag = true;
    end
    
    ctrl_cond = contains( {Conditions.name}, "Control puff", "IgnoreCase", true );
    ptx_cond = contains( {Conditions.name}, 'musc', 'IgnoreCase', true );
    if sum(ptx_cond)==0
        fprintf(1, 'Did not find PTX condition!! Continuing!!')
        fprintf(1, '%s', exp_path )
        continue
    end
    consCond = find( or(ctrl_cond, ptx_cond) );
    Nccond = length( consCond );
    prmSubs = nchoosek( 1:Nccond, 2 );

    pairedStimFlags = arrayfun(@(c) any( ...
        Conditions(1).Triggers(:,1) == ...
        reshape( Conditions(c).Triggers(:,1), 1, [] ), 2 ), consCond, fnOpts{:} );
    pairedStimFlags = cat( 2, pairedStimFlags{:} );

    consCondNames = string( { Conditions( consCond ).name  } );

    [behRes, behFig_path, behData, aInfo] = analyseBehaviour( beh_path, ...
        "ConditionsNames", cellstr( consCondNames ), ...
        "PairedFlags", pairedStimFlags, ...
        "FigureDirectory", figure_path, ...
        "ResponseWindow", [25, 350] * m, ...
        "ViewingWindow", [-450, 500] * m, ...
        "figOverWrite", fowFlag );

    if ~exist( "fr", "var" ) && ldFlag
        try
            load( get_full_path( dir( fullfile( beh_path, "RollerSpeed*.mat" ) ) ), "fr")
        catch
            load( fullfile( beh_path, 'RollerFrameRate.mat' ), 'fr' )
        end
        ldFlag = false;
    end

    [pAreas, ~, behAreaFig] = createBehaviourIndex(behRes);
    behMeasures = string({behAreaFig.Name});
    biFigPttrn = behMeasures+"%s";
    biFigPttrn = arrayfun(@(s) sprintf(s, sprintf(" %s (%%.3f)", ...
        consCondNames ) ), biFigPttrn );

    for it = 1:numel(behMeasures)
        behRes = arrayfun(@(bs, ba) setfield( bs, ...
            strrep( behMeasures(it), " ", "_" ), ba), behRes(:), pAreas(:,it) );
    end

    arrayfun(@(f) set( f, 'UserData', behRes ), behAreaFig );

    biFN = arrayfun(@(s) sprintf( biFigPttrn(s), pAreas(:,s) ), 1:numel(behMeasures) );

    arrayfun(@(f, fn) saveFigure(f, fullfile(behFig_path, fn), true, fowFlag), ...
        behAreaFig(:), biFN(:) );

    close all
end

%% 
% batch_dir = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch13_beh" ); 
% ss_fp = dir( fullfile( batch_dir, "PTX\WT*\*PTX\Behaviour\BehaviourResults*.mat" ) );
batch_dir = fullfile( "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch12_ephys.e" );
ss_fp = dir( fullfile( batch_dir, 'MC', 'vGlut*', '*PTX*', '*', 'BehaviourResults V-0.45 - 0.50 s R25.00 - 350.00 ms.mat' ) );
Nf = numel( ss_fp );
bl_changes = zeros( Nf, 8 );
my_zscore = @(x, m, s) ( x - m ) ./ ( s .* (s~=0) + 1 .* (s==0) );
for cf = 1:Nf
    load( get_full_path( ss_fp(cf) ), "behRes" )
    c_idx = contains( {behRes.ConditionName}, 'control puff', 'IgnoreCase', true );
    p_idx = contains( {behRes.ConditionName}, 'ptx', 'IgnoreCase', true );
    c = {behRes(c_idx).Results.Baseline};
    c = cat( 1, c{:} );
    t = {behRes(p_idx).Results.Baseline};
    t = cat( 1, t{:} );
    [cz, centre, scale] = zscore( c , 0, 2 );
    bl_changes(cf,:) = median( my_zscore( t, centre, scale ), 2 ) - ...
        median( cz, 2 );
end
bs_n = string( {behRes(1).Results.BehSigName} );
