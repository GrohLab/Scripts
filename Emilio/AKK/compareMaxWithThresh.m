function [moveFlag, thCell] = compareMaxWithThresh(mx, cmx)
%COMPAREMAXWITHTHRESH compares each column of the maximum value per trial
%per signal against a given threshold set.
%   Detailed explanation goes here, later
%% 
mdl = fit_poly( [1, 64], [double(cmx)/64, double(cmx)], 1 );
thCell = { ( ( 1:64 )'.^[1,0] ) * mdl };
moveFlag = mx(:) > thCell{1}(:)';
end
