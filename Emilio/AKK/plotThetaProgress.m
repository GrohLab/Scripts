function [probFigs] = plotThetaProgress(logMat, thSet, sNames)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%%
probFigs = gobjects(size(logMat,2),1);
for cs = 1:size(logMat,2)
    probFigs(cs) = figure; plot(thSet{cs}, sum(logMat{cs})/size(logMat{cs},1),...
        "DisplayName", sNames(cs))
    lgnd = legend("show"); set(lgnd, "Box", "off", "Location", "best");
    ylabel("Trial proportion"); set(gca, "Box", "off", "Color", "none");
    ylim([0,1]);
end
end

