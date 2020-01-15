%% try dpca with jittering data

%Z:\Jesus\Jittering\Already Curated\190619_Jesus_Emilio_Jittering_3703_1500_1520

clear all
dataDir='Z:\Jesus\Jittering\Already Curated\190702_Jittering_3720_1520_1520'
cd(dataDir)
load('jesus_emilio_jittering1_PopulationSpikeAnalysis.mat')
load('jesus_emilio_jittering1_all_channels.mat')
load(['jesus_emilio_jittering1analysis.mat'],'Triggers')

%% another example with gamma
clear all
dataDir='\\lsdf02.urz.uni-heidelberg.de\sd19B001\Jesus\Jittering\Already Curated\190712_Jesus_Emilio_Jittering_3700_1520_1500'
cd(dataDir)
load('Jitter_PopulationSpikeAnalysis.mat')
load('Jitter_all_channels.mat')
load(['Jitteranalysis.mat'],'Triggers')



%% another example
clear all
dataDir='Z:\Jesus\Jittering\Already Curated\190703_Jittering_3720_1520_1520'
cd(dataDir)
%load('Jitter_PopulationSpikeAnalysis.mat')
load('H5_1.3mW_all_channels.mat')
load(['H5_1.3mWanalysis.mat'],'Triggers', 'Conditions')
%% make it simple and find everything 
% well isolated clusters, only isolated
goodsInd = cellfun(@(x) x==1,sortedData(:,3));multiInd=[];
%multiInd = cellfun(@(x) x==2,sortedData(:,3));
goods = [find(goodsInd); find(multiInd)];
Ns = min(structfun(@numel,Triggers));
% Total duration of the recording
Nt = Ns/fs;  %seconds
clusters=sortedData(goods,2);ids=sortedData(goods,1);
%clusters = cellfun(@(x) round(x*fs/1000),clusters,'UniformOutput',false);  %spikes now at ms resolution
ppms=1;
biSpikes=spalloc(floor(Nt*1000*ppms),numel(goods),numel(cell2mat(clusters)));% preallocate sparse matrix to save space
% ms x num clusters
T=size(biSpikes,1); %recording length in ms



stimCond=Conditions([3:9]);
%stimCond(end-1).Triggers=stimCond(end-1).Triggers(1:5:end-1,:)
%make big matrix of all spike times
for j=1:numel(goods)
    tt=round(clusters{j}*1000*ppms);tt=tt(tt>0 & tt<=T);%spike times in ms*ppms
    S=zeros(T,1);S(tt)=1;
    biSpikes(:,j)=sparse(S);
end
%%
maxTrials=max(cellfun(@(x) size(x,1),{stimCond.Triggers}));
%in ms because of division before, though triggered seg takes samples
timeBefore=-500*ppms;
timeAfter=500*ppms;
trialLength=-timeBefore+timeAfter+1;

%preallocate omg could be huge!
% firingRates: N x S x T x maxTrialNum
firingRates=zeros(numel(clusters),numel(stimCond),trialLength,maxTrials);
for n =1:numel(clusters)
    for s=1:numel(stimCond)
        triggers=round(stimCond(s).Triggers(:,1)/fs*1000*ppms); %triggers in ms*ppms
        X=TriggeredSegments(biSpikes(:,n),triggers,timeBefore,timeAfter);
        firingRates(n,s,1:size(X,1),1:size(X,2))=(shiftdim(X,-2));  %some fuckery to deal with 4d array
    end
  
end

%all spikes for 1 neuron and 1 condition
%temp=squeeze(firingRates(1,1,1,:,:))';

trialNum=ones(numel(clusters), numel(stimCond))*maxTrials;
firingRatesAverage=mean(firingRates,4);

time = 1:trialLength;
%

    %bin triggered spikes for each neuron (use reshape business from JP)
    %squish into firingRates
    
    % trialNum: N x S x D
    % firingRates: N x S x D x T x maxTrialNum
    % firingRatesAverage: N x S x D x T

    % N is the number of neurons
    % S is the number of stimuli conditions (F1 frequencies in Romo's task)
    % D is the number of decisions (D=2)
    % T is the number of time-points (note that all the trials should have the
    % same length in time!)
    
    %include balanced number of trials for all conditions
    %do dPCA just for lag?  test if we include w and l as well
    %let's first just do PCA?
    
   
    
% Define parameter grouping

% *** Don't change this if you don't know what you are doing! ***
% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% For two parameters (e.g. stimulus and time, but no decision), we would have
% firingRates array of [N S T E] size (one dimension less, and only the following
% possible marginalizations:
%    1 - stimulus
%    2 - time
%    [1 2] - stimulus/time interaction
% They could be grouped as follows: 

combinedParams = {{1, [1 2]}, {2}};
margNames = {'Stimulus', 'Time'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
% They are marked on the plots with vertical lines
timeEvents =-timeBefore;

% check consistency between trialNum and firingRates   %%remove decision
% loop
for n = 1:size(firingRates,1)
    for s = 1:size(firingRates,2)
        %for d = 1:size(firingRates,3)
            %assert(isempty(find(isnan(firingRates(n,s,d,:,1:trialNum(n,s,d))), 1)), 'Something is wrong!')
            assert(isempty(find(isnan(firingRates(n,s,:,1:trialNum(n,s))), 1)), 'Something is wrong!')
        %end
    end
end


% setting random number of repetitions for each neuron and condition
ifSimultaneousRecording = true;  % change this to simulate simultaneous 
                                 % recordings (they imply the same number 
                                 % of trials for each neuron)

% Step 1: PCA of the dataset
close all
X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W, ...
    'combinedParams', combinedParams);

% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default, ...
   'combinedParams', combinedParams);

% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams);
toc

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);



% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'combinedParams', combinedParams, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);



%% visualize
n=1:size(firingRates,1);
t = [1:length(time)]'; %time points of interest
X = firingRatesAverage(n,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = W(:,whichMarg==1)'*X;gi
Zfull = reshape(Z(1:3,:), [3 size(firingRatesAverage,2) length(t)]);
ZZ = mean(Zfull, 3);
ZZ = Zfull;  %I THINK??? as we are missing a dimension

%

      
i=7;
id=ones(numel(t),1);
controlw=[squeeze(ZZ(1,i,:)) squeeze(ZZ(2,i,:)) squeeze(ZZ(3,i,:)) t id*i];

for i=1:5
    str=['lag' num2str(i) '=[squeeze(ZZ(1,i,:)) squeeze(ZZ(2,i,:)) squeeze(ZZ(3,i,:)) t  id*i];'];
    eval(str);
end


src = cat(1, controlw, lag1, lag2, lag3, lag4, lag5);
%make pleasant colors
[colormap]=cbrewer('seq', 'BuPu', 5)
cmap=repmat([1 0 0], numel(t),1);
for i=1:5
    cmap=[cmap; repmat(colormap(i,:), numel(t),1)];
end

src =[src brighten(cmap,0)];
comet3n(src, 'speed', 10, 'headsize', 2, 'tailwidth', 1.5, 'taillength', 50)
%%
% ------------- 3D visualization ------------

X = firingRatesAverage(n,:,t);
X = X(:,:);
X = bsxfun(@minus, X, mean(X,2));
Z = W(:,whichMarg==1)'*X;
Zfull = reshape(Z(1:3,:), [3 size(firingRatesAverage,2) length(t)]);
colors = jet(8);
colors = colors([1 2 4 5 6 7 8 ],:);
colors(3,:) = colors(3,:)*0.6;
colors(4,:) = [238 210 2]/256;
ZZ = mean(Zfull, 3);
ZZ = Zfull;  %I THINK??? as we are missing a dimension
t1 = [41:110];
t2 = [391:460];

figure('Position', [2000 0 1300 1300])
axes('Units', 'Pixels', 'Position', [1 1 1300 1300], 'Color', [0 0 0])
hold on
for i=1:6 %for each condition
    plot3(squeeze(ZZ(1,i,:)), squeeze(ZZ(2,i,:)), squeeze(ZZ(3,i,:)), 'LineWidth', 2, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t1)), squeeze(ZZ(2,i,t1)), squeeze(ZZ(3,i,t1)), 'LineWidth', 4, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t2)), squeeze(ZZ(2,i,t2)), squeeze(ZZ(3,i,t2)), 'LineWidth', 4, 'Color', colors(i,:))
    plot3(squeeze(ZZ(1,i,t2(1:2:end))), squeeze(ZZ(2,i,t2(1:2:end))), squeeze(ZZ(3,i,t2(1:2:end))), '.', ...
        'MarkerSize', 25, 'Color', colors(i,:))
end
axis normal
axis(110*[-1 1 -1 1 -1 1])
axis vis3d


%%
t2=t;
t1=t;
figure
    grid on
for i=[6 7 1:5] %for each condition
    plot3(squeeze(ZZ(1,i,:)), squeeze(ZZ(2,i,:)),squeeze(ZZ(3,i,:)),'.-', 'LineWidth', 2, 'Color', colors(i,:))
    %plot3(squeeze(ZZ(1,i,:)), squeeze(ZZ(2,i,:)),squeeze(ZZ(3,i,:)),'.-', 'MarkerSize', 10, 'Color', colors(i,:))

    grid on
    xlabel dpca1
    ylabel dpca2
    zlabel dpca3
%set(gca,'Zlim',zlims)
%set(gca,'Xlim',xlims)
%set(gca,'Ylim',ylims)
hold on
    pause
end


%%
zlims=get(gca,'Zlim')
xlims=get(gca,'Xlim')
ylims=get(gca,'Ylim')


%%
plot(squeeze(ZZ(1,i,t1)), squeeze(ZZ(2,i,t1)), 'LineWidth', 4, 'Color', colors(i,:))
plot(squeeze(ZZ(1,i,t2)), squeeze(ZZ(2,i,t2)), 'LineWidth', 4, 'Color', colors(i,:))
plot(squeeze(ZZ(1,i,t2(1:2:end))), squeeze(ZZ(2,i,t2(1:2:end))), '.', ...
        'MarkerSize', 25, 'Color', colors(i,:))


    

    
%% attempting to make a movie!
pca1=squeeze(ZZ(1,i,:));
pca2=squeeze(ZZ(2,i,:));
X=pca1
Y=nan(numel(X),1,numel(t));
Y(:,1,:)=pca2;
[hL] = animator(pca1,pca2,'ylim','auto','title','dPCA','smooth','on');
animator(pca1,pca2)
%%
    %plot3, ,squeeze(ZZ(3,i,:)),'.-', 'LineWidth', 2, 'Color', colors(i,:))
    
i=7;
controlw=[squeeze(ZZ(1,i,:)) squeeze(ZZ(2,i,:)) squeeze(ZZ(3,i,:))]
i=1
lag1=[squeeze(ZZ(1,i,:)) squeeze(ZZ(2,i,:)) squeeze(ZZ(3,i,:))]
all={};all{1}=controlw;all{2}=lag1;
plot_coords(all);
      trajectory_plotter(['test'], 10, controlw,lag1);
      
  

      
      %%
      
%%
     
  %plot 3 5-D random walks, each in a different color.  Save the output to
  %'trajectories_test.avi' in the current working directory.
  x1 = cumsum(randn(1000, 5), 1);
  x2 = cumsum(randn(1000, 5), 1);
  x3 = cumsum(randn(1000, 5), 1);
  trajectory_plotter('trajectories_test', 10, x1, x2, x3);

 

%%
vidObj = VideoWriter('practice.avi');
open(vidObj);
thetas = 1:1:360;
for i = 1:length(thetas)
    view([cosd(-135+thetas(i)) sind(-135+thetas(i)) 0.2])
    pause(0.01)
    thisFrame = getframe(gca);
    thisFrame.cdata = thisFrame.cdata(1:1127,1:1127,:);
    writeVideo(vidObj, thisFrame);
end
close(vidObj);

figure('Position', [2000 0 1300 1300])
axes('Units', 'Pixels', 'Position', [1 1 1300 1300], 'Color', [0 0 0])
hold on
axis(110*[-1 1 -1 1])
cols = [];
for i=1:7
    cols = [cols; linspace(colors(i,1),0,30)' linspace(colors(i,2),0,30)' linspace(colors(i,3),0,30)'];
end
colormap(cols)






%% Optional: estimating "signal variance"

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams, ...
    'Cnoise', Cnoise, 'numOfTrials', trialNum);

% Note how the pie chart changes relative to the previous figure.
% That is because it is displaying percentages of (estimated) signal PSTH
% variances, not total PSTH variances. See paper for more details.

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);

%% Optional: decoding

decodingClasses = {[(1:S)' (1:S)'], repmat([1:2], [S 1]), [], [(1:S)' (S+(1:S))']};

accuracy = dpca_classificationAccuracy(firingRatesAverage, firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], [], [], decodingClasses)

accuracyShuffle = dpca_classificationShuffled(firingRates, trialNum, ...
    'lambda', optimalLambda, ...
    'combinedParams', combinedParams, ...
    'decodingClasses', decodingClasses, ...
    'simultaneous', ifSimultaneousRecording, ...
    'numRep', 5, ...        % increase to 100
    'numShuffles', 20, ...  % increase to 100 (takes a lot of time)
    'filename', 'tmp_classification_accuracy.mat');

dpca_classificationPlot(accuracy, [], accuracyShuffle, [], decodingClasses)

componentsSignif = dpca_signifComponents(accuracy, accuracyShuffle, whichMarg);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16,                ...
    'componentsSignif', componentsSignif);
