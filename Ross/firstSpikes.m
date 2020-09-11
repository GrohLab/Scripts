function  firstSpikes(relativeSpkTmsStruct, gclID, dataDir)

for a = 1: length(relativeSpkTmsStruct)
    spVals(a).means = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
    spVals(a).sd = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
    spVals(a).sZcf = zeros(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
    cf = length(relativeSpkTmsStruct(a).SpikeTimes(1,:))*ones(length(relativeSpkTmsStruct(a).SpikeTimes(:,1)),1);
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
        spVals(a).sZcf(b,1) = sum(vals == 0);
        vals(find(vals == 0)) = [];
        spVals(a).means(b,1) = mean(vals);
        spVals(a).sd(b,1) = std(vals);
        cf(b) = cf(b) - spVals(a).sZcf(b,1);
    end
    x = [1:length(relativeSpkTmsStruct(a).SpikeTimes(:,1))];
    y = log10(spVals(a).sd*1000);
    f = figure('Name',[relativeSpkTmsStruct(a).name,' Standard Deviations of First Spikes Post-TTL']);
    scatter(x,y,(cf+1)*2.5,[0,0,0])
    t = text(x,y,gclID,'FontSize', 6);
    title(relativeSpkTmsStruct(a).name);
    xlabel('Unit No.');
    ylabel('log10(Std Dev) (ms)');
    savefig(f, fullfile(dataDir,[relativeSpkTmsStruct(a).name,' Standard Deviations of First Spikes Post-TTL.fig']), 'compact');
end
end


