% IDEA: 
% One big figure where we have 4 subplots of the different conditions. with
% the different puff intensities in every plot.

%putting for different puff intensities a graph together with the relevant general
% probab 
%load manually the 3-5 different intensities and plot them together 

% WHAT 6 PUFF INTENSITIES DO YOU WANT :
% puffInt = ["0.0 bar" , "0.8 bar", "1.2 bar", "1.6 bar", "2.0 bar", "2.4bar"]; WT28
puffInt = ["0.0 bar" , "0.6 bar", "1.2 bar", "1.8 bar", "2.4 bar", "3.0 bar"]; %WT27


% CHANGE FOR THE RELEVANT ANIMAL FILE 
% fName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\"; %WT28
fName = "Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\"; %WT27

% CHANGE FOR THE RELATIVE ANIMAL (#28)
% puff1 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211209\0.0bar\BehaviourData2021-12-09T19_11_48.mat");
% puff2 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211210\0.8bar\BehaviourData2021-12-10T09_49_14.mat");
% puff3 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211210\1.2bar\BehaviourData2021-12-10T10_06_27.mat");
% puff4 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211210\1.6bar\BehaviourData2021-12-10T10_16_02.mat");
% puff5 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211213\2.0bar\BehaviourData2021-12-13T16_52_44.mat");
% puff6 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211213\2.4bar\BehaviourData2021-12-13T17_06_42.mat");

puff1 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT28\211209\0.0bar\BehaviourData2021-12-09T19_11_48.mat");
puff2 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211206\0.6bar\BehaviourData2021-12-06T18_07_54.mat");
puff3 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211206\1.2bar\BehaviourData2021-12-06T18_29_29.mat");
puff4 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211207\1.8bar\BehaviourData2021-12-07T18_47_48.mat");
puff5 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211207\2.4bar\BehaviourData2021-12-07T19_08_17.mat");
puff6 = load("Z:\Emilio\SuperiorColliculusExperiments\Roller\Batch3\WT27\211207\3.0bar\BehaviourData2021-12-07T19_26_46.mat");


%% getting the differrent arrays needed for the plot !! MANUALLY
%beautify the code: (loops etc.)

logMat = cell(4,6);

logMat{1} = puff1.behResults{5};
logMat{5} = puff2.behResults{5};
logMat{9} = puff3.behResults{5};
logMat{13} = puff4.behResults{5};
logMat{17} = puff5.behResults{5};
logMat{21} = puff6.behResults{5};

logMat{2} = puff1.behResults{6};
logMat{6} = puff2.behResults{6};
logMat{10} = puff3.behResults{6};
logMat{14} = puff4.behResults{6};
logMat{18} = puff5.behResults{6};
logMat{22} = puff6.behResults{6};

logMat{3} = puff1.behResults{7};
logMat{7} = puff2.behResults{7};
logMat{11} = puff3.behResults{7};
logMat{15} = puff4.behResults{7};
logMat{19} = puff5.behResults{7};
logMat{23} = puff6.behResults{7};

logMat{4} = puff1.behResults{8};
logMat{8} = puff2.behResults{8};
logMat{12} = puff3.behResults{8};
logMat{16} = puff4.behResults{8};
logMat{20} = puff5.behResults{8};
logMat{24} = puff6.behResults{8};

%% stays the same

thSet = puff1.thSet;
sNames = ["Ipsilateral Whiskerpatch", "Contralateral Whiskerpatch", "Nose", "Rollerspeed"];
xArray = [repmat("Angle [°]",1,3), "Speed [cm/s]"];
sigTh = [5; 5; 2; 2.5];

%% FIGURE

% One big figure where we have 4 subplots of the different conditions. with
% the different puff intensities in every plot.
newcolors = [0.0 0.0 0.5977
    0.0 0.3984 0.7969
    0 0.5 1
    0.3984 0.3984 0.9961
    0 0.7 0.7
    0.25 0.90 0.7];

pfig = figure;
colororder(newcolors);
set(gcf,'Name','Different Puffintensities for 4 conditions','NumberTitle', ...
    'off','Units','normalized','Position',[0 0 1 1]); % or 'WindowState','fullscreen'
% one fig with the different puffintensities NOT different plots
for cs = 1:size(logMat,1) % cs is for the 4 different conditions (subplot)
    subplot(2,2,cs)
    for pv = 1:size(logMat,2) %puffvalues
        plot(thSet{cs}, sum(logMat{cs,pv})./size(logMat{cs,pv},1),...
            "DisplayName", puffInt(pv), 'Linewidth', 1.5);
        hold on
    end
    [t, s] = title(sprintf('Trials crossing \\theta: %s', sNames(cs)), ...
        '\sigma:' + string(sigTh(cs)),'Fontsize',16);
    s.FontAngle = 'italic';
    s.FontSize = 14;
    xlabel(sprintf('%s', xArray(cs)), 'Fontsize',14);
    ylabel("Trial proportion", 'Fontsize', 14);
    set(gca, "Box", "off", "Color", "none");
    ylim([0,1]);
    lgnd = legend("show");
    set(lgnd, "Box", "off", "Location", "best");
end
hold off

%% Saving the Figure

figname = 'PuffintensitiesWith2.4bar';
saveFigure(pfig, fullfile(fName, figname), 1);
