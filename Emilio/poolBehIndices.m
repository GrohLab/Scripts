%#ok<*AGROW,*SAGROW> 
fnOpts = {'UniformOutput', false};
expandName = @(x) fullfile(x.folder, x.name);

batchDir = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch11_ephys.MC";

animalDirs = dir(batchDir);
pointFlag = arrayfun(@(x) any(strcmpi(x.name, {'.','..'})), animalDirs);
nonAnFlag = arrayfun(@(x) isempty(regexp(x.name, ...
    '[a-zA-Z]+[0-9]+\S','match')), animalDirs);
animalDirs(pointFlag | nonAnFlag) = [];
oldMouse = "";
mc = 0; mice = [];
for cad = animalDirs(:)'
    currMouse = cad.name;
    if string(oldMouse) ~= string(currMouse)
        oldMouse = currMouse;
        mice = [mice; struct('Name', currMouse, 'Sessions',[])];
        mc = mc + 1;
        sc = 0; oldSess = "";
    end
    sessDirs = getSubFolds(expandName(cad));
    % Just date sessions
    onlyDateSessFlag = regexp(arrayfun(@(x) string(x.name), sessDirs), ...
        '[0-9]{6}', 'match');
    sessDirs(isempty(onlyDateSessFlag)) = [];
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
        ctrFlag = contains(consCondNames, 'control','IgnoreCase',true);
        mscFlag = contains(consCondNames, 'muscimol','IgnoreCase',true);
        if string(oldSess) ~= string(currSess)
            oldSess = currSess; auxStruct = struct('Date', currSess, ...
                'Control', behRes(ctrFlag).BehIndex, ...
                'Muscimol', behRes(mscFlag).BehIndex);
            if ~isfield(mice, 'Sessions')
                mice(mc).Sessions = auxStruct;
            else
                mice(mc).Sessions = [mice(mc).Sessions; auxStruct];
            end
            sc = sc + 1;
        end
    end
end