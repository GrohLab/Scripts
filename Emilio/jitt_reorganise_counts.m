cfg_struct = configStructure;
cfg_struct.Viewing_window_s = [-0.3,0.2];
[PSTHpupt, psthTx, n_trials] = getPSTH_perU_perT( ...
    relativeSpkTmsStruct, cfg_struct);
n_conditions = numel(PSTHpupt);
n_elements = cellfun(@(x) numel(x), PSTHpupt);
condition_names = {'laser1ms','laser10ms','laser50ms',...
    'laser100ms','laser200ms','whisker'};
for ccond = 1:n_conditions
    [n_bins, n_neurons] = size(PSTHpupt{ccond},[2,3]);
    
    
end