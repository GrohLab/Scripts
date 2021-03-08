% ScatterFormat
a = 5;
fa = figure(a);
% Does the figure contain one graph? If so then
axfa = fa.Children;
axfa(2).XLim = [0,60];
axfa(2).YLim = [0, 60];
pe = axfa.Children;
peClasses = arrayfun(@class,pe,'UniformOutput',0);
textType = peClasses{5};
deleteIdx = arrayfun(@(x,classType) isa(x,classType),pe,repmat(string(textType),size(pe)));
pe(deleteIdx).delete
pe = axfa.Children;
pe(5).Marker = 'o';
pe(5).CData = [0,0,0];
pe(5).SizeData = 75;
pe(5).DisplayName = 'Cluster';
pe(4).DisplayName = 'y = x';
pe(4).LineStyle = '-';
pe(3).SizeData = 50;
pe(3).CData = [1, 0, 0];
pe(3).Marker = '+';
pe(3).DisplayName = 'Significant p<0.05';
pe(1).DisplayName = 'Shuffled p<0.05';
pe(1).Marker = 'x';
pe(1).SizeData = 50;
pe(1).CData = [1, 0, 0];
configureFigureToPDF(figure(a));

% Does the figure contain multiple graphs? Then
b = 4;
axfa = fa.Children;
pe = axfa(b).Children;
axfa(b).XLim = [0, 60];
axfa(b).YLim = [0, 60];
peClasses = arrayfun(@class,pe,'UniformOutput',0);
textType = peClasses{5};
deleteIdx = arrayfun(@(x,classType) isa(x,classType),pe,repmat(string(textType),size(pe)));
pe(deleteIdx).delete
pe = axfa(b).Children;
pe(4).Marker = 'o';
pe(4).CData = [0,0,0];
pe(4).SizeData = 75;
pe(4).DisplayName = 'Cluster';
pe(3).DisplayName = 'y = x';
pe(3).LineStyle = '-';
XDataPre = pe(3).XData(1,2);
XDataPost = axfa(b).XLim(1,2);
factor = XDataPost/XDataPre;
pe(3).XData(1,2) = XDataPre*factor;
pe(3).YData(1,2) = XDataPre*factor;
pe(2).SizeData = 50;
pe(2).CData = [1, 0, 0];
pe(2).Marker = '+';
pe(2).DisplayName = 'Significant p<0.05';
m = pe(1).YData(1,2)- pe(1).YData(1,1);
f = pe(1).XData(1,2); 
g = m/f;
pe(1).YData(1,2) = (pe(1).YData(1,1) + g*pe(3).XData(1,2));
pe(1).XData(1,2) = pe(3).XData(1,2);
pe(1).DisplayName = ['Trend line ', num2str(g)];
configureFigureToPDF(figure(a));