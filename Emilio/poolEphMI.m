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
total_var_dist = @(dmat) integral( @(x) abs( pdf( dmat(1), x ) - pdf( dmat(2), x ) ), -5, 5 );
tocol = @(x) x(:);
%% Assuming 1 level of animal organisation i.e.
% BatchX/FolderA/Animal001
% BatchX/FolderB/Animal002
batchDir = fullfile( "Z:\Emilio\SuperiorColliculusExperiments", ...
    "Roller", "Batch17_ephys.MC");
%Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch15_ephys

childFolders = dir(batchDir);

pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), childFolders);
fileFlag = ~[childFolders.isdir]';
childFolders(pointFlag | fileFlag) = [];
animalFolders = arrayfun(@(d) recursiveFolderSearch(expandName(d), ...
    rsOpts{:}), childFolders, fnOpts{:}); animalFolders = cat(1, animalFolders{:});
%% Looping animals
oldMouse = "";
mc = 0; mice = [];
for cad = animalFolders(:)'
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
        miIdxFiles = arrayfun(@(d) dir(fullfile(d, "LogMI 2-conditions*.fig")), ...
            miFigDir, fnOpts{:}); miIdxFiles = cat(1, miIdxFiles{:});
        if isempty(miIdxFiles)
            fprintf(1, 'No ephys analysis done! Skipping %s!\n', curDir)
            continue
        end
        miIdxFig = arrayfun(@(x) openfig(expandName(x), 'invisible'), ...
            miIdxFiles); 
        MI = arrayfun(@(x) get(x, 'UserData'), miIdxFig, fnOpts{:} );
        arrayfun(@close, miIdxFig); 
        switch class(MI{1})
            case 'cell'
                condNames = string( cat( 1, MI{:}(1:2) ) );
                condNames = join( [condNames(1), "vs", condNames(2)] );
                miVal = cat( 1, MI{:}{3} );
            case 'struct'
                condNames = cellfun( @(m) m.MI.Comparative, MI );
                miVal = cellfun( @(m) m.MI.Value, MI );
        end
        Nbix = numel(MI); c = 1;
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
btchName = regexp( batchDir, 'Batch\d+','match' );
behFP = fullfile( batchDir, btchName+"_EphMI.mat" );
svOpts = {'-mat'};
if exist(behFP, "file")
    svOpts = {'-append'};
end
save(behFP, "mice", svOpts{:})
habFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "multi", ...
    m.Sessions), mice, fnOpts{:});
cat( 1, habFlag{:} )
%{
%% multiple
jittDist = makedist('Normal', 'mu', 0, 'sigma', 1/9);
habFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "multi", ...
    m.Sessions), mice, fnOpts{:});
habTable = arrayfun(@(m, f) {m.Sessions(f{:}).DataTable}, mice, habFlag, ...
    fnOpts{:});
pBehIdx = cellfun(@(x) cellfun(@(y) cell2mat(y.BehaviourIndices), x, ...
    fnOpts{:}), habTable, fnOpts{:});
Ncc = cellfun(@(x) cellfun(@(y) numel(y), x), pBehIdx, fnOpts{:});
rSz = cellfun(@(x) max(cellfun(@(y) numel(y), x)), pBehIdx);
cSz = cellfun(@numel, pBehIdx);
resBehIdx = arrayfun(@(x,y) nan(x,y), rSz, cSz, fnOpts{:});
resTable = cell(numel(mice), 1); Nm = numel(mice);
mNames = arrayfun(@(m) m.Name, mice); clrMap = roma(Nm);
habFig = figure('Name', 'Intensity v.s. index', 'Color', 'w');
ax = axes('Parent', habFig, 'Color', 'none', 'Box', 'off', 'NextPlot', 'add');
x = []; y = [];
for m = 1:Nm
    mxSub = find(Ncc{m} == rSz(m), 1, "first");
    for ci = 1:cSz(m)
        endS = numel(pBehIdx{m}{ci});
        resBehIdx{m}(1:endS,ci) = pBehIdx{m}{ci};
    end
    resTable{m} = table(resBehIdx{m}, ...
        'RowNames', mice(m).Sessions(mxSub).DataTable.Row, ...
        'VariableNames', "BehaviourIndices");
    x = [x; reshape(ones(rSz(m), cSz(m)).*(1:rSz(m))', [], 1)];
    y = [y; resTable{m}.BehaviourIndices(:)];
    %scatter(ax, ones(cSz(m),rSz(m)).*(1:rSz(m)) + ...
    scatter(ax, (1:rSz(m)) + random(jittDist, [1,rSz(m)]), ...
        mean(resTable{m}.BehaviourIndices,2,'omitnan')', [], clrMap(m,:), ...
        "filled", "MarkerFaceAlpha", 0.5)
end
xticks(ax, 1:max(rSz));

lgObj = legend(ax, mNames);
set(lgObj, "Box", 'off', 'Color', 'none', 'Location', 'best', 'AutoUpdate', 'off')

%% single
singFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), mice, fnOpts{:});
behTable = arrayfun(@(m, f) {m.Sessions(f{:}).DataTable}, mice, singFlag, ...
    fnOpts{:}); behTable = arrayfun(@(t) cat(1, t{:}{:}), behTable, fnOpts{:});
behTable = cat(1, behTable{:});
ctrl = behTable{behTable.Conditions == "Control Puff", "BehaviourIndices"};
ptx = behTable{behTable.Conditions == "PTX", "BehaviourIndices"};
figure; scatter(ones(size(ptx, 1),2).*[1,2], [ctrl, ptx])
hold on; plot(ones(2,size(ptx, 1)).*[1;2], [ctrl, ptx]', 'k:')
behTable = [behTable; mice(5).Sessions(2).DataTable]
ctrl = behTable{behTable.Conditions == "Control Puff", "BehaviourIndices"}
ptx = behTable{behTable.Conditions == "PTX", "BehaviourIndices"}
figure; scatter(ones(size(ptx, 1),2).*[1,2], [ctrl, ptx])
hold on; plot(ones(2,size(ptx, 1)).*[1;2], [ctrl, ptx]', 'k:')
ptx = behTable{contains(behTable.Conditions, "PTX"), "BehaviourIndices"}
figure; scatter(ones(size(ptx, 1),2).*[1,2], [ctrl, ptx])
xlim([0,3])
xticks(1:2)
hold on; plot(ones(2,size(ptx, 1)).*[1;2], [ctrl, ptx]', 'k:')
[p, h] = ranksum(ctrl, ptx)
[p, h] = ranksum(ctrl(setdiff(1:6,3)), ptx(setdiff(1:6,3)))
koFlag = true(size(ctrl));
koFlag(3) = false;
[p, h] = ranksum(ctrl(koFlag), ptx(koFlag))
[ctrl, ptx]
[ctrl, ptx, koFlag]
[p, h] = ranksum(ctrl(koFlag), ptx(koFlag), "tail", "right")
[p, h] = ranksum(ctrl(koFlag), ptx(koFlag), "tail", "left")
[p, h] = ranksum(ctrl, ptx, "tail", "left")
[p, h] = ranksum(ctrl(koFlag), ptx(koFlag), "tail", "left")
ylim([0,1])
ylabel('Behaviour index')
xticks(1:2)
xticklabels({'Control', 'PTX'})
hold on; plot([1,2], max([ctrl, ptx], [], "all")*([1,1]+0.1), 'k')
text(1.5, max([ctrl, ptx],[], "all")*1.1, '\ast', "HorizontalAlignment", 'center', "VerticalAlignment", "bottom")
title(["PTX [60 \muM] in SC";"Significance: left tail"])
configureFigureToPDF(gcf)
figure; scatter(ones(sum(koFlag),2).*[1,2], [ctrl(koFlag), ptx(koFlag)])
hold on; plot([1,2], max([ctrl(koFlag), ptx(koFlag)], [], "all")*([1,1]+0.1), 'k')
hold on; plot(ones(2,sum(koFlag)).*[1;2], [ctrl(koFlag), ptx(koFlag)]', 'k:')
xlim([0,3])
xticks(1:2)
xticklabels({'Control', 'PTX'})
ylim([0,1])
ylabel('Behaviour index')
title(["PTX [60 \muM] in SC";"Significance: left tail"])
configureFigureToPDF(gcf)
saveFigure(gcf, fullfile("Z:\Emilio\SuperiorColliculusExperiments\Roller\GenFigures", "PTX effect"), true);
text(1.5, max([ctrl, ptx],[], "all")*1.1, '\ast', "HorizontalAlignment", 'center', "VerticalAlignment", "bottom")
%%
muscFlag = arrayfun(@(m) arrayfun(@(s) cellfun(@(c) ...
    any(contains(c, 'musc', ctOpts{:}),2), s.DataTable.Conditions), ...
    m.Sessions, fnOpts{:}), mice, fnOpts{:});
sessFlag = cellfun(@(f) cellfun(@any, f), muscFlag, fnOpts{:});
behTable2 = arrayfun(@(m, f1) m.Sessions(f1{:}).DataTable, ...
    mice, sessFlag, fnOpts{:});
multFlag = cellfun(@(t) ~isstring(t.Conditions), behTable2);
behTableM = cellfun(@(t, f, s) t(f{s},:), behTable2(multFlag), ...
    muscFlag(multFlag), sessFlag(multFlag), fnOpts{:});
behTableM = cellfun(@(t) table(t.Conditions{:}(:), t.BehaviourIndices{:}(:), ...
    'VariableNames', t.Properties.VariableNames), behTableM, fnOpts{:});
behTable2 = cat(1, behTableM{:}, behTable2{~multFlag});

dateFlag = arrayfun(@(m) arrayfun(@(s) ~contains(fieldnames(s), 'Date'), ...
    m.Sessions, fnOpts{:}), mice, fnOpts{:});
mCatg = arrayfun(@(mn) categorical(regexp(mn.Name, '[A-Za-z]{2}', ...
    'match')), mice); 
Nfn = cellfun(@(x) cellfun(@sum, x), dateFlag, fnOpts{:});
values = arrayfun(@(m) arrayfun(@(s) struct2cell(s), m.Sessions, ...
    fnOpts{:}), mice, fnOpts{:});
%}