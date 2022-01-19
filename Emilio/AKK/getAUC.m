function [genProb] = getAUC(logMat)
%GETAUC does something that we will describe
%   Detailed explanation goes here, later
genProb = cellfun(@(x) nnz(x)/numel(x), logMat);
end