%#ok<*AGROW,*SAGROW> 
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);
animalPattern = '[A-Za-z]+\d{2,}';
rsOpts = {animalPattern, 'SearchType', 'expression'};
ctOpts = {'IgnoreCase', true};
lsOpts = {'L\d+.\d+', 'match'};
%% Assuming 1 level of animal organisation i.e.
% BatchX/FolderA/Animal001
% BatchX/FolderB/Animal002
batchDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC";
childFolders = dir(batchDir);

pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), childFolders);
childFolders(pointFlag) = [];
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
        ephDir = expandName(dir(fullfile(expandName(csd), 'ephys*')));
        if isempty(ephDir)
            fprintf(1, "This session didn't have ephys!\n")
            fprintf(1, "Skipping\n")
            continue
        end
        currSess = char(regexp(csd.name, '\d{6}', 'match'));
        behFigDir = recursiveFolderSearch(ephDir,...
            "Beh V-0.25 - 0.50 s R5.00 - 400.00 ms");
        behIdxFiles = dir(fullfile(behFigDir, "BehIndex*.fig"));
        behIdxFig = arrayfun(@(x) openfig(expandName(x), 'invisible'), ...
                behIdxFiles);
        behRes = arrayfun(@(x) get(x, 'UserData'), behIdxFig, fnOpts{:});
        brSz = cellfun(@numel, behRes); c = 1;
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
        if string(oldSess) ~= string(currSess)
            oldSess = currSess; auxStruct = struct('Date', currSess, ...
                'Control', behRes(ctrFlag).BehIndex, ...
                'Laser_continuous', behRes(xor(delFlag,frqFlag)).BehIndex);
            if any(frqFlag)
                auxStruct.Laser_freq = behRes(frqFlag).BehIndex;
            end
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
    end
end