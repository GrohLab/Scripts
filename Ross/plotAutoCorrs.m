% run Neurocorrs on units first!
tmReach = 300;
ms = round(length(corrs{1,1})/tmReach);
test = zeros(1,tmReach);
test(1) = sum(corrs{1,1}(1,1:ms));
for i = 2:tmReach
test(i) = sum(corrs{1,1}(1,(i-1)*ms+1:i*ms));
end
figure; plot(test)
ax = gca;
ax.XTick = [0, 50, 100, 125, 130, 135, 140, 145, 149, 150, 151, 155, 160, 165, 170, 175, 200, 250, 300];
ax.XTickLabel = ax.XTick - 150;