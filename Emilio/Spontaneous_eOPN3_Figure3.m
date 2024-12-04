%%
fnOpts = {'UniformOutput', false};
my_xor = @(x) xor( x(:,1), x(:,2) );
ovwFlag = false;

getMI = @(x,d) diff(x, 1, d)./sum(x, d).*(sum(x,d)>0) + (1.*(sum(x,d)==0 | sum(x,d)< 1e-12));
% getMI = @(ctr, cnd) (cnd - ctr) ./ ( (cnd + ctr).*((cnd + ctr)>0) + (1.*((cnd + ctr)==0 | (cnd + ctr)<1e-9) ));
load(['C:\Users\jefe_\seadrive_root\Emilio U\Für meine Gruppen\GDrive ' ...
    'GrohLab\Projects\00 SC\SC Behaviour\Figures\Figure 3\Matlab ' ...
    'figures\Data\Spontaneous FR.mat'])
expandName = @(x) fullfile( x.folder, x.name );

sessions_per_mouse = arrayfun(@(x) dir( fullfile( expandName( x ) ) ), ...
    animalFolders, fnOpts{:} );
sessions_per_mouse = cellfun(@(x) x([x.isdir] & ...
    ~ismember( {x.name}, {'.','..'} )), sessions_per_mouse, fnOpts{:} );
Nsess = sum( cellfun(@numel, sessions_per_mouse ) );
% 1. Mouse and session number, 2. spike count per spontaneous window, 3.
% spontaneous window limits, 4. firing frequency per window per unit, 5.
% treatment start, 6. conditions indicator, 7. median firing rate per unit,
% 8. p-value of modulation index
results = cell(Nsess, 8 );
cr = 1;
for ca = 1:numel(animalFolders)
    fprintf(1, 'Animal: %s\n', animalFolders(ca).name )
    sessions = sessions_per_mouse{ca};
    ce = 1;
    for cs = 1:numel(sessions)
        results{cr,1} = [ca,cs];
        fprintf(1, 'Session: %s\n', sessions(cs).name );
        spkFN = dir( fullfile( expandName( sessions(cs) ), 'ephys*', ...
            '*spike_times.mat' ) );
        anFN = dir( fullfile( expandName( sessions(cs) ), '*', ...
            '*analysis.mat' ) );
        if ~isempty( spkFN )
            load( expandName(spkFN ) )
        else
            fprintf(1, "%s: ", sessions(cs).name )
            fprintf(1, 'No spike_times file found!\n')
            fprintf(1, 'Skipping!\n')
            continue
        end
        if ~isempty( anFN )
            load( expandName(anFN ), 'Conditions', 'fs' )
        else
            fprintf(1, "%s: ", sessions(cs).name )
            fprintf(1, 'No analysis file found!\n')
            fprintf(1, 'Skipping!\n')
            continue
        end
        FigureDir = fullfile( spkFN.folder, "Figures" );
        Nu = numel( spike_times );
        Topn = Conditions(2).Triggers(1,1)/fs;
        Texp = Conditions(4).Triggers(end,2)/fs;
        allTrigs = sortrows( cat(1, Conditions(2:4).Triggers ), 1, "ascend" );

        lmts = [[0; allTrigs(:,2)/fs + 0.2], [allTrigs(:,1)/fs;Texp]];
        lmts(diff( lmts, 1, 2 ) < 0,:) = [];

        treat_sub = find( lmts(:,1) < Topn, 1, "last" );
        Nt = size(lmts, 1);
        Nsp = cellfun(@(u) arrayfun(@(t) sum( my_xor( u < lmts(t,:) ) ), ...
            1:size( lmts, 1 ) ), spike_times, fnOpts{:} );
        Nsp = cat( 1, Nsp{:} );

        fr_sp = Nsp ./ diff( lmts, 1, 2 )';

        cndID = [ones( 1, treat_sub), 1+ones(1, Nt-treat_sub)];
        cndID = repmat( cndID, Nu, 1 );

        medFr = arrayfun(@(c) median( fr_sp(:,cndID(1,:)==c), 2 ), 1:2, fnOpts{:} );
        medFr = cat(2, medFr{:});

        medMI = getMI( medFr, 2 );
        p = signrank( medMI );

        results(cr,2:8) = {Nsp, ... Spike cound per window per unit
            lmts, ... Window limits
            fr_sp, ... Firing rate per unit per window
            treat_sub, ... Subscript when treatment started
            cndID(1,:), ... Condition membership indicator 
            medFr, ... Median firing rate per unit over each condition windows
            [p, median( medMI( ~isnan( medMI ) ) )] }; % p-value and median MI for all units per condition.
        cr = cr + 1;
        %%

        f = figure("Color", "w"); t = createtiles( f, 1, 4 );
        ax(1) = nexttile(t,[1,3]);

        loglog(ax(1), medFr(:,1), medFr(:,2), "k." );
        hold(ax(1), 'on'); line(ax(1), xlim, xlim, 'Color', 'k', 'LineStyle', ':' )
        xticklabels(ax(1), xticks(ax(1)) ); yticklabels(ax(1), yticks(ax(1)) )
        xlabel(ax(1), 'Control [Hz]', 'interpreter', 'latex' );
        ylabel(ax(1), 'eOPN3 [Hz]', 'interpreter', 'latex'  )
        axis( ax(1), 'square' )
        ytickangle(ax(1), 90 )

        ax(2) = nexttile(t);
        boxchart(ax(2), getMI(medFr, 2), 'Notch', 'on', ...
            'BoxFaceColor', 0.15*ones(1,3), 'JitterOutliers', 'on', ...
            'MarkerStyle', '.', 'MarkerColor', 0.15*ones(1,3) )
        text( ax(2), 1, -1.1, sprintf("$p=%.3g$", p), 'Interpreter', 'latex', ...
            "HorizontalAlignment", "center" )
        ylabel(ax(2), 'Decrease $\leftarrow$ Modulation index $\rightarrow$ Increase', ...
            'Interpreter', 'latex' )
        yline( ax(2) , 0, 'k:' )
        disappearAxis(ax(2))
        ylim( ax(2), [-1,1] )
        cleanAxis(ax);
        set( ax, 'TickDir', 'out' )
        set( f, 'UserData', {cndID(1,:), medFr, medMI, p} )
        title(t, 'Spontaneous firing rate per unit', 'interpreter', 'latex')
        saveFigure(f, fullfile( FigureDir, 'Spontaneous firing rate pre and post eOPN3' ), true, ovwFlag )
    end
    close all
    fprintf(1, 'Complete!\n' )
end
%%
results( all( cellfun(@isempty, results ), 2 ), : ) = [];
rTable = cell2table( results, "VariableNames", {'ID', 'SpikeCounts', 'SpontaneousWindows', ...
    'FR', 'TreatmentStart', 'ConditionID', 'MedianFR_pu', 'p_MedianTot'} );
save( "c:\Users\jefe_\seadrive_root\Emilio U\Für meine Gruppen\GDrive " + ...
    "GrohLab\Projects\00 SC\SC Behaviour\Figures\Figure 3\Matlab " + ...
    "figures\Data\Spontaneous FR.mat", "rTable", "animalFolders" )