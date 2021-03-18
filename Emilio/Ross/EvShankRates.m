% EvShankRatePercent
ControlInd = find(clInfo.ActiveUnit == true);
rF = responseWindow(2) - responseWindow(1);
clInfo{ControlInd,'Rate_Ev_Control'} = mean(Counts{1,2}')'/(responseWindow(2) - responseWindow(1));
clInfo{ControlInd,'Rate_Ev_CNO'} = mean(Counts{2,2}')'/(responseWindow(2) - responseWindow(1));
e1 = find(clInfo.CNO_MR == true &  clInfo.shank == 1);
e2 = find(clInfo.CNO_MR == true & clInfo.shank == 2);
e3 = find(clInfo.CNO_MR == true & clInfo.shank == 3);
e4 = find(clInfo.CNO_MR == true & clInfo.shank == 4);
e5 = find(clInfo.CNO_MR == true & clInfo.shank == 5);
e6 = find(clInfo.CNO_MR == true & clInfo.shank == 6);
e1Control = clInfo.Rate_Ev_Control(e1);
e1CNO = clInfo.Rate_Ev_CNO(e1);
e2Control = clInfo.Rate_Ev_Control(e2);
e2CNO = clInfo.Rate_Ev_CNO(e2);
e3Control = clInfo.Rate_Ev_Control(e3);
e3CNO = clInfo.Rate_Ev_CNO(e3);
e4Control = clInfo.Rate_Ev_Control(e4);
e4CNO = clInfo.Rate_Ev_CNO(e4);
e5Control = clInfo.Rate_Ev_Control(e5);
e5CNO = clInfo.Rate_Ev_CNO(e5);
e6Control = clInfo.Rate_Ev_Control(e6);
e6CNO = clInfo.Rate_Ev_CNO(e6);
for a = 1:length(e1Control)
e1R(a,1) = ((e1CNO(a,1) - e1Control(a,1))/e1Control(a,1))*100;
end
for a = 1:length(e2Control)
e2R(a,1) = ((e2CNO(a,1) - e2Control(a,1))/e2Control(a,1))*100;
end
for a = 1:length(e3Control)
e3R(a,1) = ((e3CNO(a,1) - e3Control(a,1))/e3Control(a,1))*100;
end
for a = 1:length(e4Control)
e4R(a,1) = ((e4CNO(a,1) - e4Control(a,1))/e4Control(a,1))*100;
end
for a = 1:length(e5Control)
e5R(a,1) = ((e5CNO(a,1) - e5Control(a,1))/e5Control(a,1))*100;
end
for a = 1:length(e6Control)
e6R(a,1) = ((e6CNO(a,1) - e6Control(a,1))/e6Control(a,1))*100;
end
figure; subplot(1, 6, 1);
fprintf([' no. Inf = ', num2str(numel(e1R(e1R == inf))), '  \n']);
e1R(e1R == inf) = [];
bar(e1R);
subplot(1, 6, 2);
fprintf([' no. Inf = ', num2str(numel(e2R(e2R == inf))), '  \n']);
e2R(e2R == inf) = [];
bar(e2R);
subplot(1, 6, 3);
fprintf([' no. Inf = ', num2str(numel(e3R(e3R == inf))), '  \n']);
e3R(e3R == inf) = [];
bar(e3R);
subplot(1, 6, 4);
fprintf([' no. Inf = ', num2str(numel(e4R(e4R == inf))), '  \n']);
e4R(e4R == inf) = [];
bar(e4R);
subplot(1, 6, 5);
fprintf([' no. Inf = ', num2str(numel(e5R(e5R == inf))), '  \n']);
e5R(e5R == inf) = [];
bar(e5R);
subplot(1, 6, 6);
fprintf([' no. Inf = ', num2str(numel(e6R(e6R == inf))), '  \n']);
e6R(e6R == inf) = [];
bar(e6R);