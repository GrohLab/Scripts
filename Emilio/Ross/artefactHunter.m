firstSpikes{a} = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),length(relativeSpkTmsStruct(a).SpikeTimes(1,:)));
artVals(a).means = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
artVals(a).sd = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
for a = 1: length(relativeSpkTmsStruct)
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
        artVals(a).means(b,1) = mean(vals);
        artVals(a).sd(b,1) = std(vals);
    end
end

