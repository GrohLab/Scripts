function [moveFlag] = compareMaxWithThresh(mx, thCell)
%COMPAREMAXWITHTHRESH compares each column of the maximum value per trial
%per signal against a given threshold set.
%   Detailed explanation goes here, later
%% 
% Validation section for size and compatibility.
% Number of signals
Ns = size(mx, 2);

% Number of given cells
Nth = numel(thCell);

if Ns ~= Nth
    fprintf(1, "The signals in the stack is different from the number of")
    fprintf(1, " sets given!\n")
    fprintf(1, "Please, check your variables!\n")
    moveFlag = {};
    return
end

moveFlag = arrayfun(@(x) mx(:,x) > thCell{x}, 1:Ns,...
    "UniformOutput", false);

end
