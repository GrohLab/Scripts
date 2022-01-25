function [moveFlag] = compareMaxWithThresh(mx, thCell)
%COMPAREMAXWITHTHRESH compares each column of the maximum value per trial
%per signal against a given threshold set.
%   Detailed explanation goes here, later
%% 
moveFlag = mx' > thCell{1};

end
