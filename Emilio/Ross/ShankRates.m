% ShankSepsRate

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
s1R(a,1) = ((e1CNO(a,1) - e1Control(a,1))/e1Control(a,1))*100;
end
for a = 1:length(e2Control)
s2R(a,1) = ((e2CNO(a,1) - e2Control(a,1))/e2Control(a,1))*100;
end
for a = 1:length(e3Control)
s3R(a,1) = ((e3CNO(a,1) - e3Control(a,1))/e3Control(a,1))*100;
end
for a = 1:length(e4Control)
s4R(a,1) = ((e4CNO(a,1) - e4Control(a,1))/e4Control(a,1))*100;
end
for a = 1:length(e5Control)
s5R(a,1) = ((e5CNO(a,1) - e5Control(a,1))/e5Control(a,1))*100;
end
for a = 1:length(e6Control)
s6R(a,1) = ((e6CNO(a,1) - e6Control(a,1))/e6Control(a,1))*100;
end
figure; subplot(1, 6, 1);
fprintf([' no. Inf = ', num2str(numel(percentdif(s1R == inf))), '  \n']);
s1R(s1R == inf) = [];
bar(s1R);
subplot(1, 6, 2);
fprintf([' no. Inf = ', num2str(numel(percentdif(s2R == inf))), '  \n']);
s2R(s2R == inf) = [];
bar(s2R);
subplot(1, 6, 3);
fprintf([' no. Inf = ', num2str(numel(percentdif(s3R == inf))), '  \n']);
s3R(s3R == inf) = [];
bar(s3R);
subplot(1, 6, 4);
fprintf([' no. Inf = ', num2str(numel(percentdif(s4R == inf))), '  \n']);
s4R(s4R == inf) = [];
bar(s4R);
subplot(1, 6, 5);
fprintf([' no. Inf = ', num2str(numel(percentdif(s5R == inf))), '  \n']);
s5R(s5R == inf) = [];
bar(s5R);
subplot(1, 6, 6);
fprintf([' no. Inf = ', num2str(numel(percentdif(s6R == inf))), '  \n']);
s6R(s6R == inf) = [];
bar(s6R);