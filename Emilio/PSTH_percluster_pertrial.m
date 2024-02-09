histOpts2 = {'BinLimits', [0,0.3], 'BinWidth', 1e-3, 'Normalization', ...
    'probability'};
binSizes = logspace(-4,-2,100);
b_var = zeros(ceil(diff(histOpts2{2})/histOpts2{4}),length(binSizes),'single'); 
ii = 1;
for bz = binSizes
    binSize = bz; vw = configStructure.Viewing_window_s;
    histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
        'Normalization', 'count'};
    fnOpts = {'UniformOutput', false};
    psth_u = arrayfun(@(v) cellfun(@(u) histcounts(u, histOpts{:}), ...
        relativeSpkTmsStruct(1).SpikeTimes(v,:), fnOpts{:}), ...
        1:size(relativeSpkTmsStruct(1).SpikeTimes,1), fnOpts{:});
    psth_u = cellfun(@(x) cat(1, x{:}), psth_u, fnOpts{:});
    psth_u = cat(3, psth_u{:});
    
    %mdl_psth_tx = fit_poly([1,size(psth_u, 2)], vw + [1,-1]*binSize/2, 1);
    %psth_tx = ( ( 1:size(psth_u, 2) )'.^[1,0] ) * mdl_psth_tx;
    %plot(psth_tx, mean(psth_u, [1,3]))
    b_var(:,ii) = histcounts(mean(psth_u, [1,3]), histOpts2{:});
    ii = ii + 1;
    % figure; imagesc( psth_tx, [], psth_u)
end
figure; imagesc(1:length(binSizes), histOpts2{2}, log10(b_var+1e-6))
axis xy; colormap(inferno)