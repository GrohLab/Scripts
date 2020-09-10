function [spVals] = firstSpikes(relativeSpkTmsStruct)


for a = 1: length(relativeSpkTmsStruct)
    spVals(a).means = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
    spVals(a).sd = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
    for b = 1: length(relativeSpkTmsStruct(a).SpikeTimes(:,1))
        vals = zeros(1,length(relativeSpkTmsStruct(a).SpikeTimes(1,:)));
        for c = 1: length(relativeSpkTmsStruct(a).SpikeTimes(1,:))
            ind = find(relativeSpkTmsStruct(a).SpikeTimes{b,c}(:,:) >= 0);
            if sum(relativeSpkTmsStruct(a).SpikeTimes{b,c} >= 0)~= false
                first = ind(1);
                vals(1,c) = relativeSpkTmsStruct(a).SpikeTimes{b,c}(first);
            else
            end
        end
        vals(find(vals == 0)) = [];
        %spVals(a).Vals(b,1) = vals;
        spVals(a).means(b,1) = mean(vals);
        spVals(a).sd(b,1) = std(vals);
    end
    figure('Name',[relativeSpkTmsStruct(a).name,' Standard Deviations of First Spikes Post-TTL']);
    plot(spVals(a).sd, 'Color',[0 0 0],'MarkerSize',10,'Marker','.',...
    'LineStyle','none');
    title(relativeSpkTmsStruct(a).name);
    xlabel('Unit No.');
    ylabel('Std Dev (s)');
end
end


