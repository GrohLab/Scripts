function [fig, TaggedIDs] = plotTaggedFraction(latencyCutoffs, sdCutoffs, clInfo, clIDs, ConditionName, TriggerTimes, sortedData, samplingFrequency)
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
Latencies = TriggerLatencies(sortedData(spkInd,2), TriggerTimes, fs, 5e-2);
mn = (cellfun(@mean, Latencies)*1e3);
sd = (cellfun(@std, Latencies)*1e3);
Dpth = -1*Dpth;
nan = isnan(mn);

latencyCutoffs = sort(latencyCutoffs, 'ascend');
sdCutoffs = sort(sdCutoffs, 'ascend');

tagged = latencyCutoffs(1) <= mn & mn <= latencyCutoffs(2) & sdCutoffs(1) <= sd & sd <= sdCutoffs(2);
taggedInd = spkInd(tagged);
TaggedIDs = sortedData(taggedInd,1);

nontagged = ~tagged | nan;

minDepth = Dpth(min(find(tagged)));
maxDepth = Dpth(max(find(tagged)));
fig = figure('Name', name, 'Color', 'White');
x = [sum(tagged), sum(~tagged)];
explode = [1, 0];
% labels = {['Tagged Fraction = ',num2str(tg), '%'], ' '};
pie(x, explode)




ax = gca;
ax.Children(2).FaceColor = [0.5 0.5 0.5];
ax.Children(4).FaceColor = [0 1 1];
ax.FontSize = 30;
% ax.FontName = 'Times New Roman';

lgd = legend;
lgd.String{1} = ['Opto-tagged'];
text(1.25,1.0,['Selection Criteria:'],'FontSize', 15 );
text(1.25,0.75,['Mean Latencies between ' num2str(latencyCutoffs(1)), ' and ',num2str(latencyCutoffs(2)), 'ms'], 'FontSize', 10);
text(1.25,0.50,['Standard Deviations between ' num2str(sdCutoffs(1)), ' and ',num2str(sdCutoffs(2)), 'ms'], 'FontSize', 10);
text(1.25,-0.25,'Opto-tagged Units Displayed:', 'FontSize', 15);
text(1.25,-0.5,['Depths between ' num2str(minDepth), ' and ',num2str(maxDepth), '\mum'], 'Interpreter', 'tex', 'FontSize', 10);
lgd.String{2} = ['Untagged '];
lgd.Location = 'northoutside';
%lgd.Box = 'off';


end


