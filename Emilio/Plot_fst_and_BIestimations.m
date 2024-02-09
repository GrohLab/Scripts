jitDist = makedist("Normal", "mu", 0, "sigma", 0.1);
vw = [-10, 50]*1e-3; binSize = 5e-4;
histOpts = {'BinLimits', vw, 'BinWidth', binSize, ...
    'Normalization', 'probability'};
uH = arrayfun(@(c) ...
    arrayfun(@(u) ...
    histcounts( (params.alpha(:,u,c) + params.g(:,c))*fs_scale + fs_centre, ...
    histOpts{:} ), 1:Ncl, fnOpts{:}), 1:Ncond, fnOpts{:});
uH = cellfun(@(c) cat(1, c{:}), uH, fnOpts{:});
uH = cat(3, uH{:});

[Ncl, Nt, Ncond] = size(uH);

mdl = fit_poly([1,Nt], vw + [1,-1]*(binSize/2), 1);
tx = ( (1:Nt)'.^[1,0] ) * mdl;

clrMap = [zeros(1,3); inferno(5)];
for cu = 1:Ncl
    figure;
    contour(ones(Nt,1)*(1:6), tx * ones(1,Ncond), squeeze(uH(cu,:,:)), 5);
    colormap(inferno(10))
    hold on;
    plot(1:6, ...
        squeeze(mean( params.alpha(:,cu,:) + reshape(params.g,2e3,1,6), 1 ) ) * ...
        fs_scale + fs_centre, "k", "LineWidth", 2)

    arrayfun(@(c) scatter( c+random( jitDist, ...
        size([firstSpkStruct(c).FirstSpikeTimes{cu,:}]) ), ...
        ...
        [firstSpkStruct(c).FirstSpikeTimes{cu,:}], ...
        "MarkerEdgeColor","k", "MarkerFaceColor",clrMap(c,:), ...
        "MarkerFaceAlpha",0.5), 1:Ncond)
    ylim([0,50]*1e-3); set(gca, "Box", "off", "Color", "none");
    yticklabels(yticks*1e3); ylabel('Latency [ms]')
    xticks(1:Ncond); xticklabels(string({relativeSpkTmsStruct.name}))
    set(get(gca, "YAxis"), "Scale", "log")
end