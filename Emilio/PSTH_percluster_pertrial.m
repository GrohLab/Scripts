binSize = 1e-4;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'count'};
fnOpts = {'UniformOutput', false};
psth_u = arrayfun(@(v) cellfun(@(u) histcounts(u, histOpts{:}), ...
    relativeSpkTmsStruct(1).SpikeTimes(v,:), fnOpts{:}), ...
    1:size(relativeSpkTmsStruct(1).SpikeTimes,1), fnOpts{:});
psth_u = cellfun(@(x) cat(1, x{:}), psth_u, fnOpts{:});
psth_u = cat(3, psth_u{:});
% psth_u = cat(1, psth_u{:});
mdl_psth_tx = fit_poly([1,size(psth_u, 2)], vw + [1,-1]*binSize/2, 1);
psth_tx = ( ( 1:size(psth_u, 2) )'.^[1,0] ) * mdl_psth_tx;
% figure; imagesc( psth_tx, [], psth_u)