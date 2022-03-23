function [moveFlag] = compareMaxWithThresh(mx, thCell)
%COMPAREMAXWITHTHRESH compares each column of the maximum value per trial
%per signal against a given threshold set.
%   Detailed explanation goes here, later
%% 
if isrow(mx)
    mx = mx';
end
moveFlag = mx > thCell{1};
end
