
Import_Pipeline %this call all of the importing code, do only once, for everything

%% pick out laser stimuli only for experiments with laser stimulation, do only once
%analysis file to save relevant parameters
%find laser triggers
analysisFile = [path '\analysis\' matfilename 'analysis']; %just the path and name as string 
load([path '\analysis\' matfilename '_denoise'],'laserSignal'); %fetch the laser signal from datas
%find laser triggers
stObj = StepWaveform(laserSignal,fs,'On/Off','Laser');
lt = stObj.Triggers; %light on/off
ltOn = find(lt(:,1));ltOff=find(lt(:,2));
fprintf(1,'Got the condition triggers\n')


laserPowers = [0, 0.15, 2, 6, 10, 14, 20, 23];
%laserPowers = [23];
% The laserTrials contains the number of pulses per laser intensity
% (Clustered into one figure)
laserTrials = 3;

close all
figure
plot(laserSignal)
hold on
plot(ltOn,laserSignal(ltOn),'ro')
save(analysisFile,'ltOn','ltOff','laserPowers','laserTrials', '-append')


%%  plot some examples
ephyData=dataCell; 

%online filtering if you want to look at particular freq. band
% Alpha: 7  - 13 Hz
% Beta:  13 - 30 Hz 
% Gamma: 30 - 70 Hz 
% Delta: 1  - 4  Hz 
% Theta: 4  - 7  Hz

%for cch = 1:size(ephyData,1)
   % fprintf(1,'Dealing with channel %d\n', cch)
   % [alpha, beta, gamma, delta, theta] = brainwaves(ephyData,ephyData{cch},fs)
   % ephyData{cch}= alpha;
%end


%UNCOMMENT THIS IF YOU WANT TO LOOK AT SPIKES
%ephyData=dataCell_sp; 


%% PARAMETERS INITIALIZATION -- Modify this section for different experiments if necessary
% The timeLapse variable contains the time before the pulse in milliseconds
% (1 seconds = 1000 milliseconds & 1 millisecond = 0.001 second) and the
% time after the pulse in milliseconds. (30 ms before and 30 ms after)
%can we automate this to figure out desired units?
m = 1e-3; %millisecond scaling factor
k = 1e3; %kilo scaling
timeLapse = [500,1500]*m;
% The laserPowers variable contains the series of laser intensities used to
% stimulate.

[~,cSt] =...
    getStacks(false(1,cols), ltOn, 'on', timeLapse, fs ,fs ,[],ephyData);

% Saving plots? if this is 'false' the figures will NOT be saved. If it is
% 'true', the figures will be saved.
saveFlag = false;
Nch = rows;
% Plotting the continuous stack
close all
timeOffset=timeLapse(1);
f = gobjects(8);
for laserPower = 1:numel(laserPowers)
    f(laserPower) =...
        figure('Name',sprintf('Power %.2f mW',laserPowers(laserPower)),...
        'Color',[1,1,1],'Renderer','painters','RendererMode','manual');
    for laserTrial = 1:laserTrials
        plotEEGchannels(...
            reshape(cSt(:,:,laserTrial + laserTrials*(laserPower - 1)),...
            Nch,sum(timeLapse)*fs+1),1:Nch,sum(timeLapse)+(1/fs),fs,0.5,timeOffset,f(laserPower));
        ax = f(laserPower).Children;
        ax.NextPlot = 'add';
    end
    
     mSignals=reshape(mean(cSt(:,:,laserPower*3 - 2:laserPower*3),3),...
        Nch,sum(timeLapse)*fs+1);
    
    [~,meanP] = plotEEGchannels(mSignals,1:32,sum(timeLapse)+(1/fs),fs,0.5,timeOffset,f(laserPower));
    set(meanP,'Color',[0,0,0])
    if saveFlag
        % figName = fullfile(path, sprintf('LaserStimulation_%.2fmW',laserPowers(laserPower))) %#ok<*UNRCH>
        figName = ['LaserStimulation_' num2str(laserPowers(laserPower))]
        savefig(f(laserPower),[figName,'.fig'])
        print(f(laserPower),[figName,'.eps'],'-depsc','-bestfit')
    end
    
    % mean triggered signal across channels, overlaid
    % this eventually should be put elsewhere--maybe some flag in plot
    % EEGchannels.
    figure
    dt = 1/fs;
    Ns = size(mSignals,2);
    timeS = ((0:Ns-1) * dt) - timeOffset;
    plot(timeS/m,mSignals);
    axis tight
    grid on
    xlabel ms!!!!
end

%%


%if triggers are available, use those
%if not, use randomly selected points