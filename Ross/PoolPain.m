% PoolPain is a script that takes in all the experiments to be analysed
% together and will then analyse them ala DE_Jittering/PAIN scripts


%%
% Select the experiment folders and add to the Obj
Obj.Experiments = struct;
nExpts = str2double(cell2mat(inputdlg('How Many Experiments are you considering?')));

for exp = 1:nExpts
    Obj.Experiments(exp).Directory = uigetdir('Z:\',...
        'Choose a working directory');
    if Obj.Experiments(exp).Directory == 0
        return
    end
end

%% Get data from the selected directories

ObjNames = {'bin', 'smrx', 'analysis', 'all_channels', 'cluster_info'};
for exp = 1:nExpts
    
    cDir = Obj.Experiments(exp).Directory;
    % fprintf(['Reading the relevant files in ', cDir , '. \n This may take a few mins'])
    files = {dir(fullfile(cDir, '*.bin')),...
        dir(fullfile(cDir, '*.smrx')),...
        dir(fullfile(cDir, '*analysis.mat')),...
        dir(fullfile(cDir, '*all_channels.mat')),...
        dir(fullfile(cDir, '*cluster_info.tsv'))};
    
    for name = 1:length(ObjNames)
        Obj.Experiments(exp).(ObjNames{name}) = {};
    end
    
    
    for ind = 1:length(files)
        file = files{ind};
        
        if length(file) > 1
            promptStr = {'Select the relevant file: '};
            input = struct2cell(file);
            input = (cell(input(1,1:end)))';
            fileNo = listdlg('PromptString', promptStr, 'ListString', input);
            file = dir(fullfile(cDir, file(fileNo).name));
            Obj.Experiments(exp).(ObjNames{ind}) = file;
        elseif length(file) < 1
            fprintf(['There is no ', ObjNames{ind}, ' file in directory ', num2str(exp), ' \n']);
            fprintf('Maybe you can find one manually? \n \n');
        else
            Obj.Experiments(exp).(ObjNames{ind}) = file;
        end
        
    end
    
    for mat = [3,4]
        indx = Obj.Experiments(exp).(ObjNames{mat});
        Obj.Experiments(exp).(ObjNames{mat}) = load([indx.folder, '\', indx.name]);
    end
    
    tsv = 5;
    indx = Obj.Experiments(exp).(ObjNames{tsv});
    Obj.Experiments(exp).(ObjNames{tsv}) = getClusterInfo([indx.folder, '\', indx.name]);
end

