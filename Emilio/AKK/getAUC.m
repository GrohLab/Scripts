function [genProb] = getAUC(logMat)
%GETAUC does something that we will describe
%   Detailed explanation goes here, later
genProb = nnz(logMat)/numel(logMat);
end