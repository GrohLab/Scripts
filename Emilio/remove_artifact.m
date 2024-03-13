trig_samples = 5;
Nts = 2*trig_samples + 1;
trig_offset_win = (-trig_samples:trig_samples);
random_triggers = sort( randsample(length(laser_triggers), samples) );
for cchan = 1:64
    for crt = 1:length(laser_triggers)
        artif_win = laser_triggers(crt) + trig_offset_win;
        segm = double( data( cchan, artif_win ) );
        correction = segm(1) * ((Nts:-1:1)/(Nts+1)) + ...
            segm(end) * ((1:Nts)/(Nts+1)) + ...
            rand(1, Nts)*5 - 2.5;

        artif_diff = (correction - segm).^2;
        corr_fin = (1 - artif_diff/max(artif_diff)).^3;
        segm2 = segm.*corr_fin + (1-corr_fin).*correction;
        segm = int16(round(segm2));
        data(cchan, artif_win) = segm;
    end
end