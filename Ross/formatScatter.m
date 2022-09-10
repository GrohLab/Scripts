function fig = formatScatter(fig)
set(fig, 'Position', [0.2220, 0.1781, 0.5446, 0.6943])
ax = findobj(fig.Children, 'Type', 'axes');
for ca = 1: length(ax)
    ax(ca).FontName = 'Arial';
    ax(ca).FontSize = 15;
    sc = findobj(ax(ca).Children, 'Type', 'Scatter');
    for csc = 1:length(sc)
        sc(csc).SizeData = 400;
%         if csc =~ 3
%         sc(csc).CData = [0, 0, 1];
%         end
    end
end
end


