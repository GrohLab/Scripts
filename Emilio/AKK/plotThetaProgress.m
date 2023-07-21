function [probFigs] = plotThetaProgress(logMat, thSet, sNames, varargin)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%
p = inputParser;

addRequired(p, 'logMat', @(x) islogical(x{:}))
addRequired(p, 'thSet', @(x) isvector(x{:}))
addRequired(p, 'sNames', @(x) ~isempty(x))
addParameter(p, 'showPlots', true, @(x) islogical(x) & numel(x) == 1)

parse(p, logMat, thSet, sNames, varargin{:})

logMat = p.Results.logMat;
thSet = p.Results.thSet;
sNames = p.Results.sNames;
showPlots = p.Results.showPlots;

probFigs = gobjects(size(logMat,2),1);
for cs = 1:size(logMat,2)
    probFigs(cs) = figure('Visible', showPlots); 
    plot(thSet{cs}, sum(logMat{cs})/size(logMat{cs},1),...
        "DisplayName", sNames(cs))
    lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best");
    ylabel("Trial proportion"); set(gca, "Box", "off", "Color", "none");
    ylim([0,1]);
end
end

