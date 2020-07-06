%% PARAMETERS INITIALIZATION -- Modify this section for different experiments if necessary
% The timeLapse variable contains the time before the pulse in milliseconds
% (1 seconds = 1000 milliseconds & 1 millisecond = 0.001 second) and the
% time after the pulse in milliseconds. (30 ms before and 30 ms after)
m = 1e-3;
k = 1e3;

timeLapse = [30, 30]
%timeLapse = [30, 30]*m; % The m means that the values are now in
%milliseconds
%timeLapse = [0.01, 0.02];
% The laserPowers variable contains the series of laser intensities used to
% stimulate.
laserPowers = [0, 10, 50, 70];
%laserPowers = [70];
% The laserTrials contains the number of pulses per laser intensity
% (Clustered into one figure)
%laserTrials = 3;
laserTrials = 100;
% Saving plots? if this is 'false' the figures will NOT be saved. If it is
% 'true', the figures will be saved. 
saveFlag = true;
% figurePath need to be changed to the folder were I want to save my
% images!
figurePath = uigetdir('F:\Kilosorting\','Where do you want to save your figures?');
if figurePath == 0
    fprintf(1,'Cancelling...\n')
    return
end

% Filtering question - which frequency band would you like to see?
continueLoop = true;
while continueLoop
    freqAns = inputdlg('Please provide a frequency band (in Hz) for band-pass filtering?',...
        'Band-pass', [1, 60], {'600 9000'});
    if isempty(freqAns)
        fprintf(1,'Cancelling...\n')
        continueLoop = false;
        return
    end
    try
        freqBand = cell2mat(cellfun(@str2num, freqAns, 'UniformOutput', false));
    catch
        fprintf(1,'Sorry, there''s an error: %s', freqAns)
    end
    if numel(freqBand) ~=2
        fprintf(1,'Sorry, input not correct')
    else
        continueLoop = false;
    end
end
fprintf(1, 'The cut frequency are the following: %.3f and %.3f\n', freqBand(1),...
    freqBand(2))
%% Loading the recording
read_Intan_RHD2000_file
fprintf(1,'Now filtering your data...\n')
fs = frequency_parameters.amplifier_sample_rate;
[rows, cols] = size(amplifier_data);
if cols < rows
     amplifier_data = amplifier_data';
     aux = rows;
     rows = cols;
     cols = aux;
end
dataCell = cell(rows,1);
for cch = 1:min(size(amplifier_data))
    fprintf(1,'Dealing with channel %d\n', cch)
    dataCell(cch) = {iir50NotchFilter(amplifier_data(cch,:),fs)};
    dataCell(cch) = {iirSpikeFilter(dataCell{cch},fs, freqBand)};
end
digChanInd = sum(board_dig_in_data,2);
plv = digChanInd > 0;
names = num2str(find(plv));
% names = mat2cell(names,[1,1]);
if sum(plv) > 1
    figure;h = plot(board_dig_in_data(plv,:)');
    legend(names)
    pause
    [choseVar, iOk] = listdlg('ListString',names,...
        'SelectionMode','single','PromptString','Choose the laser variable');
    if iOk
        LaserSub = str2double(names(choseVar));
    else
        fprintf(1,'Cancelling...')
        return
    end
else
    LaserSub = find(plv);
end
fprintf(1,'The digital channel for the laser will be %d\n',LaserSub)
laserSignal = board_dig_in_data(LaserSub,:);
clearvars -except laser* fs timeLapse saveFlag cols rows m k dataCell figurePath
Nch = rows;
fprintf(1,'Finished loading data\n')
%% Getting the condition trigger
% These variables detect the rising and falling edges of the digital input
% pulse (laser).
stObj = StepWaveform(laserSignal,fs,'On/Off','Laser');
clearvars laserSignal
% The lt variable contains a logical array size 2(on/off) x
% total_number_of_samples indicating the beginning and end of a pulse with
% a true value.
lt = stObj.Triggers;
% This contains the sample number at which the pulse rose. If you want to
% change for pulse fall, write as follows: ltOn = find(lt(:,2));
ltOn = find(lt(:,1));
clearvars lt
[~,cSt] =...
    getStacks(false(1,cols), ltOn, 'on', timeLapse * m, fs ,fs ,[],dataCell);
fprintf(1,'Got the condition triggers\n')

%% Plotting the continuous stack 
f = gobjects(8);
for laserPower = 1:numel(laserPowers)
    f(laserPower) =...
        figure('Name',sprintf('Power %.2f mW',laserPowers(laserPower)),...
        'Color',[1,1,1],'Renderer','painters','RendererMode','manual');
    for laserTrial = 1:laserTrials
        plotEEGchannels(...
            reshape(cSt(:,:,laserTrial + laserTrials*(laserPower - 1)),...
            Nch,sum(timeLapse*1e-3)*fs+1),1:Nch,sum(timeLapse*1e-3)+(1/fs),...
            fs,0.5,timeLapse(1)*m,f(laserPower));
        ax = f(laserPower).Children;
        ax.NextPlot = 'add';
    end
    [~,meanP] = plotEEGchannels(...
        reshape(mean(cSt(:,:,laserPower*3 - 2:laserPower*3),3),...
        Nch,sum(timeLapse*1e-3)*fs+1),1:32,sum(timeLapse*1e-3)+(1/fs),fs,...
        0.5,timeLapse(1)*m,f(laserPower));
    set(meanP,'Color',[0,0,0])
    if saveFlag
        figName = fullfile(figurePath,...
            sprintf('LaserStimulation_%.2fmW',laserPowers(laserPower))); %#ok<*UNRCH>
        savefig(f(laserPower),[figName,'.fig']) 
        print(f(laserPower),[figName,'.eps'],'-depsc')
    end
end



%% 

