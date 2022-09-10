name = 'Corrected_H5_ChanMap';
load('Z:\Ross\ProbeFiles\Adaptor-Amplifier_Mapping.mat')


% x = ones(64,1);
% x(2:2:end) = 1+22.5;
% y = linspace(12.5,800,64)';
% depths = flip([0:20:1275]');

intan_ordered_channelmap = NaN(64,1);
% intan_ordered_depths = NaN(64,1);
for pin = 0:63
    ind = find(rhdInput == pin);
    intan_ordered_channelmap(pin+1) = electrodes_depth_ordered_H3(ind);
%     intan_ordered_depths(pin+1) = depths(ind);
end
%%
chMap = NaN(64,3);
for ind = 1:64
    ch = intan_ordered_channelmap(ind);
    dpth = find(h5chans_x_y(:,1)==ch);
    chMap(ind,:) = h5chans_x_y(dpth,:);
end

% chanMap = [1:64]';
% chanMap0ind = chanMap-1;
% xcoords = ones(64,1);
% ycoords = intan_ordered_depths;
% kcoords = ones(64,1);
% connected = true(64,1);
% fs = 1/33.3e-6;
% save('Z:\Ross\ProbeFiles\Corrected_H5_ChanMap.mat','chanMap','chanMap0ind','xcoords','ycoords','kcoords','name','connected', 'fs');