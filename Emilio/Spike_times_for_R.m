fnOpts = {'UniformOutput', false};
% N_spikes_pcond = arrayfun(@(c) sum(arrayfun(@(u) ...
%     numel([relativeSpkTmsStruct(c).SpikeTimes{u,:}]), ...
%     1:size(relativeSpkTmsStruct(c).SpikeTimes,1))), ...
%     1:length(relativeSpkTmsStruct));
N_spikes_ptrial_u25 = arrayfun(@(c) cellfun(@(ut) numel(ut), ...
    relativeSpkTmsStruct(c).SpikeTimes(25,:)), ...
    1:length(relativeSpkTmsStruct), fnOpts{:});
% (Neuron), condition, trial, spike time. No neuron information now because
% trying only with unit X (25).
unit_spike_times = zeros(sum([N_spikes_ptrial_u25{:}]), 3);

dat_i = 1; dat_j = 1;
for ccond = 1:length(relativeSpkTmsStruct)
    for cu = 25
        cond_id = repmat(ccond, sum(N_spikes_ptrial_u25{ccond}), 1);
        for ctr = 1:size(relativeSpkTmsStruct(ccond).SpikeTimes,2)
            if N_spikes_ptrial_u25{ccond}(ctr)
                trial_id = repmat(ctr, N_spikes_ptrial_u25{ccond}(ctr), 1);
                curr_i = dat_i:N_spikes_ptrial_u25{ccond}(ctr)+dat_i-1;
                unit_spike_times(curr_i,2:3) = [trial_id, ...
                    relativeSpkTmsStruct(ccond).SpikeTimes{cu,ctr}'];
                dat_i = dat_i + N_spikes_ptrial_u25{ccond}(ctr);
            else
                continue
            end
        end
        unit_spike_times(dat_j:sum(N_spikes_ptrial_u25{ccond})+dat_j-1,1) = ...
            cond_id; dat_j = dat_j + sum(N_spikes_ptrial_u25{ccond});
    end
end