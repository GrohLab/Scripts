for a = 1:length(Conditions)
    AllLaser(a,1) = contains(Conditions(a).name, 'All_Laser', 'IgnoreCase', true);
end
if sum(AllLaser) ~=1
    fprintf('Outputs for All_Laser do not = 1....Check Conditions! \n')
    return
end
ind = find(ismember(sortedData(:,1),ID))';
laserTriggers = Conditions(AllLaser).Triggers(:,1)/fs;
FirstSpikes = cell((length(ind)),length(laserTriggers));
r = 1;
for a = ind
    for b = 1:length(laserTriggers) - 1
        spikes = sortedData{a,2} >= laserTriggers(b) & sortedData{a,2} < laserTriggers(b+1);
        if sum(spikes) == 0
            FirstSpikes{r,b} = [];
        else
            spikes = find(spikes);
            first = spikes(1);
            FirstSpikes{r,b} = sortedData{a,2}(first) - laserTriggers(b);
        end
    end
    
    spikes = sortedData{a,2} >= laserTriggers(end);
    if sum(spikes) == 0
        FirstSpikes{r,b} = [];
    else
        spikes = find(spikes);
        first = spikes(1);
        FirstSpikes{r,b} = sortedData{a,2}(first) - laserTriggers(end);
    end
    r = r + 1;
end
MeanFirsts = NaN(length(ind),1);
SDcells = NaN(length(ind),1);
for a = 1:length(ind)
    FSrow = cat(2,FirstSpikes(a,1:end));
    MeanFirsts(a) = mean(cell2mat(FSrow));
    SDcells(a) = std(cell2mat(FSrow));
end
x = [1:length(SDcells)];
y = log10(MeanFirsts*1000);
scale = SDcells*10;
figure('Name', 'FirstSpikes', 'Color', 'White')
scatter(x, y, scale, [0,0,0])