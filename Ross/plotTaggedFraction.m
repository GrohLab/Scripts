function fig = plotTaggedFraction(latencyCutoffs, sdCutoffs, depthCutoffs, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)
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
name = ConditionName;
name(strfind(name, '_')) = ' ';
name = ['Tagged Fraction: ', name];
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;

latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');
% depthCutoffs =[];
if isempty(depthCutoffs)
    depthCutoffs = [min(Dpth), max(Dpth)];
elseif depthCutoffs > 0
    depthCutoffs = -depthCutoffs;
end
depthCutoffs = sort(depthCutoffs, 'ascend');


tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2)...
    & depthCutoffs(1) <= Dpth & Dpth <= depthCutoffs(2);

fig = figure('Name', name, 'Color', 'White');
x = [sum(tagged), length(mn)-sum(tagged)];
explode = [1, 0];
% labels = {['Tagged Fraction = ',num2str(tg), '%'], ' '};
pie(x, explode)




ax = gca;
ax.Children(2).FaceColor = [0 0 0];
ax.Children(4).FaceColor = [0 1 1];
ax.FontSize = 30;
% ax.FontName = 'Times New Roman';

lgd = legend;
lgd.String{1} = ['Opto-tagged'];
text(1.25,0.75,['Mean Latency between ' num2str(latencyCutoffs(1)), ' and ',num2str(latencyCutoffs(2)), 'ms']);
text(1.25,0.50,['Standard Deviation between ' num2str(sdCutoffs(1)), ' and ',num2str(sdCutoffs(2)), 'ms']);
text(1.25,0.25,['Depth between ' num2str(depthCutoffs(2)), ' and ',num2str(depthCutoffs(1)), '\mum'], 'Interpreter', 'tex');
lgd.String{2} = ['Untagged '];
lgd.Location = 'northoutside';

end


