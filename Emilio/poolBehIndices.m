%#ok<*AGROW,*SAGROW> 
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);
animalPattern = '[A-Za-z]+\d{1,}';
rsOpts = {animalPattern, 'SearchType', 'expression'};
ctOpts = {'IgnoreCase', true};
lsOpts = {'L\d+.\d+', 'match'};
behFF = "Beh V-0.25 - 0.50 s R5.00 - 400.00 ms";
tblOpts = {'VariableNames', {'Conditions', 'BehaviourIndices'}};
%% Assuming 1 level of animal organisation i.e.
% BatchX/FolderA/Animal001
% BatchX/FolderB/Animal002
batchDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC";
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
    [~, currMouse] = fileparts(cad);
    if string(oldMouse) ~= string(currMouse)
        oldMouse = currMouse;
        mice = [mice; struct('Name', currMouse, 'Sessions',[])];
        mc = mc + 1;
        sc = 0; oldSess = "";
    end
    sessDirs = getSubFolds(cad);
    % Just date sessions
    onlyDateSessFlag = arrayfun(@(x) string(regexp(x.name, '[0-9]{6}', ...
        'match')), sessDirs, fnOpts{:}); 
    sessDirs(cellfun(@isempty, onlyDateSessFlag)) = [];
    for csd = sessDirs(:)'
        curDir = expandName(csd);
        currSess = char(regexp(csd.name, '\d{6}', 'match'));
        childFolders = getSubFolds(curDir);
        sessOrgDirs = arrayfun(@(d) string(d.name), childFolders);
        sessOrgDirs(contains(sessOrgDirs, {'behaviour', 'ephys', 'figures'}, ...
            ctOpts{:})) = [];
        behFigDir = arrayfun(@(d) recursiveFolderSearch(expandName(d), ...
            behFF), childFolders, fnOpts{:}); behFigDir = cat(1, behFigDir{:});
        %{
        if isempty(ephDir)
            fprintf(1, "This session didn't have ephys!\n")
            fprintf(1, "Skipping\n")
            continue
        else
            behFigDir = recursiveFolderSearch(ephDir,...
                "Beh V-0.25 - 0.50 s R5.00 - 400.00 ms");
        end
        ephDir = expandName(ephDir);
        %}
        behIdxFiles = arrayfun(@(d) dir(fullfile(d, "BehIndex*.fig")), behFigDir);
        behIdxFig = arrayfun(@(x) openfig(expandName(x), 'invisible'), ...
                behIdxFiles); behRes = arrayfun(@(x) get(x, 'UserData'), ...
                behIdxFig, fnOpts{:}); arrayfun(@close, behIdxFig)
        condNames = cellfun(@(x) arrayfun(@(y) string(y.ConditionName), x), ...
            behRes, fnOpts{:}); behIdx = cellfun(@(x) arrayfun(@(y) ...
            y.BehIndex, x), behRes, fnOpts{:});
        brSz = cellfun(@numel, behRes); c = 1; Nbix = numel(brSz);
        sessType = 'single';
        if Nbix == 1
            dataTable = table([condNames{:}]', [behIdx{:}]', tblOpts{:});
        elseif Nbix > 1
            % We need to check where are all of these different
            % measurements are coming from.
            sessType = 'multi';
            if numel(sessOrgDirs) == Nbix
                % Same folders and measurements. Ideal situation for
                % several measurements.
                dataTable = table(condNames, behIdx, tblOpts{:}, ...
                    'RowNames', sessOrgDirs);
            else
                dataTable = table(condNames, behIdx, tblOpts{:});
            end
        end
        if string(oldSess) ~= string(currSess)
            oldSess = currSess; 
            auxStruct = struct('Date', currSess, ...
                'DataTable', dataTable, 'Type', sessType);
            %{
            if any(xor(delFlag, frqFlag))
                auxStruct.Laser_continuous = behRes(xor(delFlag, ...
                    frqFlag)).BehIndex;
            end
            if any(frqFlag)
                auxStruct.Laser_freq = behRes(frqFlag).BehIndex;
            end
            if any(musFlag)
                auxStruct.Muscimol = behRes(musFlag).BehIndex;
            end
            if any(ptxFlag)
                auxStruct.PTX = behRes(ptxFlag).BehIndex;
            end
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            %}
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
        %{
        if numel(behIdxFiles) > 1
            fprintf(1, "Found more than 1 file! Estimating the correct\n")
            Ncond = arrayfun(@(m) arrayfun(@(s) numel(fieldnames(s))-1, ...
                m.Sessions), mice, fnOpts{:}); Ncond = cat(1, Ncond{:});
            [~, whr] = min(Ncond - brSz',[],"all"); 
            if std(brSz)
                [~, c] = ind2sub([numel(Ncond), numel(brSz)], whr);
            else
                posCCN = cellfun(@(bc) arrayfun(@(bs) ...
                    string(bs.ConditionName), bc), behRes, fnOpts{:});
                c = cellfun(@(nms) any(contains(nms, 'delay', ctOpts{:})), posCCN);
            end
            fprintf(1, "Chose the following conditions:\n")
            fprintf(1, " - %s\n", arrayfun(@(x) string(x.ConditionName), ...
                behRes{c}))
        elseif isempty(behIdxFig)
            fprintf(1, "Found no BehIndex figure!\n")
            continue
        end
        behRes = behRes{c};
        arrayfun(@close, behIdxFig)
        consCondNames = arrayfun(@(x) string(x.ConditionName), behRes);
        % TODO: fix the condition names and run the loop
        ctrFlag = contains(consCondNames, 'control', ctOpts{:});
        delFlag = contains(consCondNames, 'delay', ctOpts{:});
        frqFlag = ~cellfun(@isempty,regexp(consCondNames, lsOpts{:}));
        musFlag = contains(consCondNames, 'muscimol', ctOpts{:});
        ptxFlag = contains(consCondNames, 'ptx', ctOpts{:});
        if string(oldSess) ~= string(currSess)
            oldSess = currSess; auxStruct = struct('Date', currSess, ...
                'Control', behRes(ctrFlag).BehIndex);
            if any(xor(delFlag, frqFlag))
                auxStruct.Laser_continuous = behRes(xor(delFlag, ...
                    frqFlag)).BehIndex;
            end
            if any(frqFlag)
                auxStruct.Laser_freq = behRes(frqFlag).BehIndex;
            end
            if any(musFlag)
                auxStruct.Muscimol = behRes(musFlag).BehIndex;
            end
            if any(ptxFlag)
                auxStruct.PTX = behRes(ptxFlag).BehIndex;
            end
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
        %}
    end
end
%%
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

%%
%{
dateFlag = arrayfun(@(m) arrayfun(@(s) ~contains(fieldnames(s), 'Date'), ...
    m.Sessions, fnOpts{:}), mice, fnOpts{:});
mCatg = arrayfun(@(mn) categorical(regexp(mn.Name, '[A-Za-z]{2}', ...
    'match')), mice); 
Nfn = cellfun(@(x) cellfun(@sum, x), dateFlag, fnOpts{:});
values = arrayfun(@(m) arrayfun(@(s) struct2cell(s), m.Sessions, ...
    fnOpts{:}), mice, fnOpts{:});
%}