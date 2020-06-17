% ShankSepsRate

clInfo{ControlInd,'Rate_Sp_Control'} = mean(Counts{1,1}')'/(responseWindow(2) - responseWindow(1));
clInfo{ControlInd,'Rate_Sp_CNO'} = mean(Counts{2,1}')'/(responseWindow(2) - responseWindow(1));
s1 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 1);
s2 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 2);
s3 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 3);
s4 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 4);
s5 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 5);
s6 = find(clInfo.Spont_CNO_R == true & clInfo.shank == 6);
s1Control = clInfo.Rate_Sp_Control(s1);
s1CNO = clInfo.Rate_Sp_CNO(s1);
s2Control = clInfo.Rate_Sp_Control(s2);
s2CNO = clInfo.Rate_Sp_CNO(s2);
s3Control = clInfo.Rate_Sp_Control(s3);
s3CNO = clInfo.Rate_Sp_CNO(s3);
s4Control = clInfo.Rate_Sp_Control(s4);
s4CNO = clInfo.Rate_Sp_CNO(s4);
s5Control = clInfo.Rate_Sp_Control(s5);
s5CNO = clInfo.Rate_Sp_CNO(s5);
s6Control = clInfo.Rate_Sp_Control(s6);
s6CNO = clInfo.Rate_Sp_CNO(s6);
for a = 1:length(s1R)
s1R(a,1) = ((s1CNO(a,1) - s1Control(a,1))/s1Control(a,1))*100;
end
for a = 1:length(s2Control)
s2R(a,1) = ((s2CNO(a,1) - s2Control(a,1))/s2Control(a,1))*100;
end
for a = 1:length(s3Control)
s3R(a,1) = ((s3CNO(a,1) - s3Control(a,1))/s3Control(a,1))*100;
end
for a = 1:length(s4Control)
s4R(a,1) = ((s4CNO(a,1) - s4Control(a,1))/s4Control(a,1))*100;
end
for a = 1:length(s5Control)
s5R(a,1) = ((s5CNO(a,1) - s5Control(a,1))/s5Control(a,1))*100;
end
for a = 1:length(s6Control)
s6R(a,1) = ((s6CNO(a,1) - s6Control(a,1))/s6Control(a,1))*100;
end
figure; subplot(1, 6, 1);
s1R(s1R == inf) = [];
bar(s1R);
subplot(1, 6, 2);
s2R(s2R == inf) = [];
bar(s2R);
subplot(1, 6, 3);
s3R(s3R == inf) = [];
bar(s3R);
subplot(1, 6, 4);
s4R(s4R == inf) = [];
bar(s4R);
subplot(1, 6, 5);
s5R(s5R == inf) = [];
bar(s5R);
subplot(1, 6, 6);
s6R(s6R == inf) = [];
bar(s6R);