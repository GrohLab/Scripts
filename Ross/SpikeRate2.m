
% This script gives spike rate bar graphs for user-specified bin width for
% all units in a recording, WHEN THE INPUT DATA IS THE SPIKE TIMES.

% We need the bin size.
promptStrings = {'Bin size [s]:'};
defaultInputs = {'0.05',};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defaultInputs);
binSz = str2double(answ(1));

   
% We also need the sampling frequency, which can change from recording
% to recording.
promptStrings = {'Sampling Frequency (fs)'};
defaultInputs = {'3.003003003003003e+04',};
answ = inputdlg(promptStrings,'Inputs', [1, 30],defaultInputs);
fs = str2double(answ(1));


% Instead of calculating a bin window and summing the samples for each window into the bin to get the counts,
% as you would for binary data, we need to get number of elements whose
% values fall between a bin size in seconds.


% Need to load CondSig file as well as sortedData to get the npoints.

% binSamples is the number of elements in spiketrain per bin.
binSamples = fs*binSz;

% Need the number of bins we're gonna pop our data into.
nBins = round(head68.npoints/binSamples);


% Need to create an empty matrix that's eagerly awaiting all our cluster
% counts.
badsIdx = cellfun(@(x) x==3,sortedData(:,3));
nBads = sum(badsIdx);
% Need to know the number of bad clusters that we not going to bother with.
szT = length(sortedData)- nBads;
counts = zeros(nBins, szT);

for b = 1:szT   % cluster by cluster
    if sortedData{b,3} == 3 % don't want to waste time on bad clusters
    else
        for a = 1:nBins
            fsV = fs*sortedData{b,2}'; 

            logicalfsV = (fsV >(a-1)*binSamples) & (fsV <= a*binSamples);
            counts(a,b) = sum(logicalfsV);
        end
        % rate = counts/binSz;
        % figure; bar(rate(:,b));
    end
end
 rate = counts/binSz;
% Rate will be given in Counts per Second.
% If you had 20 counts in a half-second bin size, the rate would be 40
% counts = 20/ bin size.
% So you're dividing counts by the bin size.

% We've churned out the firing rate per cluster fo a given bin width.
% Lets plot them next to eachother to get population activity.
% Experiment time in msecs for the RateMap x-axis.
exptT = [0:nBins]'*(binSz/0.001);
 
figure; imagesc(exptT, [1:szT], rate');
% x-axis = exptT, y-axis = [1:szT], plot = rate'.