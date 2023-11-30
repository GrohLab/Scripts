getMI = @(a, c) (c - a) ./ (a + c);
getRC = @(a, c) (c - a) ./ c;
normDist = makedist("Normal", "mu", 0, "sigma", 0.075);

singFlag = arrayfun(@(m) arrayfun(@(s) string(s.Type) == "single", ...
    m.Sessions), mice, fnOpts{:});
behTable = arrayfun(@(m, f) {m.Sessions(f{:}).DataTable}, mice, singFlag, ...
    fnOpts{:});

behMI = cellfun(@(m) cellfun(@(t) arrayfun(@(c) ...
    getMI(t{1,"BehaviourIndices"}, t{c,"BehaviourIndices"}), ...
    2:size(t,1)), ...
    m, fnOpts{:}), ...
    behTable, fnOpts{:});

behMI = cellfun(@(m) cat(1, m{:}), behMI, fnOpts{:});

behRC = cellfun(@(m) cellfun(@(t) arrayfun(@(c) ...
    getRC(t{1,"BehaviourIndices"}, t{c,"BehaviourIndices"}), ...
    2:size(t,1)), ...
    m, fnOpts{:}), ...
    behTable, fnOpts{:});

behRC = cellfun(@(m) cat(1, m{:}), behRC, fnOpts{:});
%%
figure; boxplot([cat(1, behMI{1:2}); padarray(behMI{3}, [0,2], nan, "both")], ...
    "Notch", "on"); hold on
title('Modulation index on behaviour index for MC\rightarrowSC excitation')
for cm = 1:3
    auxBMI = padarray(behMI{cm}, [0,(6 - size(behMI{cm},2))/2], nan, "both");
    xpos = reshape(repmat(1:size(auxBMI,2), size(auxBMI, 1 ) , 1 ) + ...
        random(normDist, size( auxBMI ) ), [], 1);
    scatter(xpos, auxBMI(:), "filled")
end
yline(0, 'k:')
xticklabels(condNames{1}{1}(2:end))
ylabel('MI(BI)'); xlabel('Conditions')
set(gca, 'Box', 'off', 'Color', 'none')
ylim([-1,1])
%%
figure; boxplot(cat(1, behMI{4:6}), "Notch", "on"); hold on
title('Modulation index on behaviour index for MC\rightarrowSC inhibition')
for cm = 4:6
    scatter(ones(size(behMI{cm}))+random(normDist, size(behMI{cm})), ...
        behMI{cm}, "filled")
end
xticklabels(condNames{4}{1}(2))
yline(0, 'k:'); set(gca, 'Box', 'off', 'Color', 'none')
