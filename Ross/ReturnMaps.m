pisi =SalSpontISI;
currentIsi = cellfun(@(x) x(1:end-1)/fs, pisi, 'UniformOutput', 0);
nextIsi = cellfun(@(x) x(2:end)/fs, pisi, 'UniformOutput', 0);
isiReturnMap = {cat(1, currentIsi{:}), cat(1, nextIsi{:})};
isiReturnMap = cellfun(@log10, isiReturnMap, 'UniformOutput', 0);
isiReturnMap = cat(2, isiReturnMap{:});
figure('Color', 'White');
Hist = histogram2(isiReturnMap(:,1), isiReturnMap(:,2), 125, 'Normalization', 'probability');
title('ISI Return Map')
xlabel('ISI_n _(_m_s_e_c_s_)');
ylabel('ISI_n_+_1 _(_m_s_e_c_s_)')
legend([{'Saline'}, {'CFA'}], 'Location', 'Northeast')
xlim([-3, 4]);
ylim([-3, 4]);
zlim([0, 0.015])
zlabel('Probability')
fig = gcf;
ax = gca;
ax.FontSize = 20;
xticks(-3:1:4)
ax.XTickLabel = (10.^cellfun(@str2double,ax.XTickLabel)) * 1e3;
ax.XTickLabelRotation = -45;
yticks(-3:1:4);
ax.YTickLabel = (10.^cellfun(@str2double,ax.YTickLabel)) * 1e3;
ax.YTickLabelRotation = 45;
