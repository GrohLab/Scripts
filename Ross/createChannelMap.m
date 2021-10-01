name = 'Corrected_H3_ChanMap';
load('Z:\Ross\ProbeFiles\Adaptor-Amplifier_Mapping.mat')

depths = flip([0:20:1275]');

intan_ordered_channelmap = NaN(64,1);
intan_ordered_depths = NaN(64,1);
for pin = 0:63
    ind = find(rhdInput == pin);
    intan_ordered_channelmap(pin+1) = electrodes_depth_ordered_H3(ind);
    intan_ordered_depths(pin+1) = depths(ind);
end



chanMap = [1:64]';
chanMap0ind = chanMap-1;
xcoords = ones(64,1);
ycoords = intan_ordered_depths;
kcoords = ones(64,1);
connected = true(64,1);
fs = 3.003003003003003e+04;
save('Z:\Ross\ProbeFiles\Corrected_H3_ChanMap.mat','chanMap','chanMap0ind','xcoords','ycoords','kcoords','name','connected', 'fs');