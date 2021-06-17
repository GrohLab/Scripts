function TaggedIDs = TagID(latencyCutoffs, sdCutoffs, clInfo, clIDs, TriggerTimes, sortedData, samplingFrequency)
fs = samplingFrequency;
clInd = ismember(clInfo.id, clIDs);
depths = table(clInfo.id(clInd), clInfo.AbsDepth(clInd));
depths = sortrows(depths,'Var2','ascend');
ID = depths{:,1};
Dpth = depths.Var2;
spkInd = [];
for i = 1:length(ID)
    spkInd = [spkInd; find(ismember(sortedData(:,1), ID(i)))];
end
tm = 5e-2;
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, tm);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);

lowestLat = 3;
lowestSD = 0.1;


mn(mn<lowestLat) = NaN;
sd(sd<lowestSD) = NaN;
nan = isnan(mn) | isnan(sd);


latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');

tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2);
taggedInd = spkInd(tagged);
TaggedIDs = sortedData(taggedInd,1);

end
