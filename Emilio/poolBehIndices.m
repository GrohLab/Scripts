%#ok<*AGROW,*SAGROW> 
%% Auxiliary variables and functions
fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);
animalPattern = '[A-Za-z]+\d{2,}';
rsOpts = {animalPattern, 'SearchType', 'expression'};
ctOpts = {'IgnoreCase', true};
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
        currSess = csd.name;
        behFigDir = recursiveFolderSearch(expandName(csd),...
            "Beh V-0.25 - 0.50 s R5.00 - 400.00 ms");
        behIdxFiles = dir(fullfile(behFigDir, "BehIndex*.fig"));
        if numel(behIdxFiles) == 1
            behIdxFig = openfig(expandName(behIdxFiles), 'invisible');
        elseif numel(behIdxFiles) > 1
            fprintf(1, "Found more than 1 file! Taking the first.\n")
            behIdxFig = openfig(expandName(behIdxFiles(1)), 'invisible');
        else
            fprintf(1, "Found no BehIndex figure!\n")
            continue
        end
        behRes = get(behIdxFig, 'UserData'); close(behIdxFig)
        consCondNames = arrayfun(@(x) string(x.ConditionName), behRes);
        % TODO: fix the condition names and run the loop
        ctrFlag = contains(consCondNames, 'control', ctOpts{:});
        delFlag = contains(consCondNames, 'delay', ctOpts{:});
        if string(oldSess) ~= string(currSess)
            oldSess = currSess; auxStruct = struct('Date', currSess, ...
                'Control', behRes(ctrFlag).BehIndex, ...
                'Muscimol', behRes(delFlag).BehIndex);
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
    end
end