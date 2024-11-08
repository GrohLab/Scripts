time_scale = 2.5e-3; delta_t = 5e-4; synchTrial = [];

for ccond = 1:length(relativeSpkTmsStruct)
    for ctr = 1:size(relativeSpkTmsStruct(ccond).SpikeTimes,2)
        spks_in_trial = [relativeSpkTmsStruct(ccond).SpikeTimes{:,ctr}];
        c_end = vw(1) + time_scale; c_init = vw(1);
        while c_end <= vw(2)
            spks_in_cons = spks_in_trial(spks_in_trial > c_init & ...
                spks_in_trial < c_end);
            if ~isempty(spks_in_cons) && numel(spks_in_cons) > 1
                dm = pdist([spks_in_cons(:)], "euclidean");
                if numel(dm) > 1
                    lnorm_fit = fitdist(dm(:), "Normal");
                    synchTrial = [synchTrial; ccond, ctr, ...
                        mean([c_init, c_end]), ...
                        log(lnorm_fit.ParameterValues(1))];
                else
                    synchTrial = [synchTrial; ccond, ctr, ...
                        mean([c_init, c_end]), log(dm(:))];
                end
            end
            c_init = c_init + delta_t;
            c_end = c_end + delta_t;
        end
    end
end