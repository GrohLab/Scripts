
figure; hold on
Nu = size(relativeSpkTmsStruct(1).SpikeTimes, 1);
clrMap = hsv(Nu);
for ccond = 1:length(relativeSpkTmsStruct)
    for cu = 15
        for ct = 1:size(relativeSpkTmsStruct(ccond).SpikeTimes, 2)
            spks_cu_ct = [relativeSpkTmsStruct(ccond).SpikeTimes{cu, ct}]';
            tid = repmat(ct, length(spks_cu_ct), 1);
            cid = repmat(cu, length(spks_cu_ct), 1);
            scatter3( spks_cu_ct, tid, cid, [], ...
                'MarkerFaceColor', clrMap(cu,:), 'MarkerEdgeColor', 'none', ...
                'MarkerFaceAlpha', 3/4)
        end
    end
end

%%
clrMap = [zeros(1,3); lines(length(relativeSpkTmsStruct)-1)];
cu = 18;
figure; hold on
tc = 1;
for ccond = 1:3
    for ctr = 1:size(relativeSpkTmsStruct(ccond).SpikeTimes, 2)
        cspks = relativeSpkTmsStruct(ccond).SpikeTimes{cu,ctr};
        if ~isempty(cspks)
            line(cspks, tc, 'Marker', 'o', ...
                'MarkerEdgeColor', 'none', ...
                'MarkerFaceColor', clrMap(ccond,:))
        end
        tc = tc + 1;
    end
end
yticks(Na/2 + cumsum([0, Na(1:end-1)]));
yticklabels({relativeSpkTmsStruct.name})
ylim([0, sum(Na)+1])