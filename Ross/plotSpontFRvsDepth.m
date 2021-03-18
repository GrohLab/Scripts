function fig = plotSpontFRvsDepth(ID, clInfo, sortedData, Triggers, samplingFrequency)

clIDs = ID;
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
Dpth = -1*Dpth;
jitter = randi(20,size(Dpth));
Dpth = Dpth + jitter; % adding small random jitter for visualisation

spkSubs = cellfun(@(x) round(x.*fs), sortedData(spkInd,2),...
'UniformOutput', false);

[firingRatePerCluster, deltaTrigTimeSum, sponSpks, sponIsi] =...
    getSpontFireFreq(spkSubs, Triggers, [0, Triggers(end,end)], fs, 2);


if class(ID) == 'cell'
    figureName = ID{1};
else
    figureName = ID;
end
figureName(strfind(figureName, '_')) = ' ';
figureName = ['Spont FR vs Depth: ', figureName];
fig = figure('Name', figureName, 'Color', 'White');

scatter(firingRatePerCluster, Dpth, 500, [1, 0, 0], '*')
xlabel('Spont. Firing Rate_{(Hz)}')
xlim([0, 100])
xticks([0:10:100])
ylim([round(round(min(Dpth),-2)/2,-2)*2-200,0]);
yticks([-1600:200:0]);
ylabel('Depth_{(\mum)}', 'Interpreter', 'tex');
ax = gca;
ax.FontSize = 15;
end
