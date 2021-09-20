
load('Z:\Ross\ProbeFiles\Adaptor-Amplifier_Mapping.mat')
intan_ordered_channelmap = NaN(64,1);
depths = flip([0:20:1275]');


for pin = 0:63
    ind = find(rhdInput == pin);
    intan_ordered_channelmap(pin+1) = electrodes_depth_ordered_H3(ind);
end

intan_ordered_depths = NaN(64,1);
for pin = 0:63
    depthInd = find(rhdInput == pin);
    intan_ordered_depths(pin+1) = depths(depthInd);
end