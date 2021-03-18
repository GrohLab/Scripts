% UnevokedByShank
ControlInd = find(clInfo.ActiveUnit == true);
rF = responseWindow(2) - responseWindow(1);
clInfo{ControlInd,'Rate_Sp_Control'} = mean(Counts{1,1}')'/rF;
clInfo{ControlInd,'Rate_Sp_CNO'} = mean(Counts{2,1}')'/rF;
s1 = find(clInfo.ActiveUnit == true & clInfo.shank == 1);
s2 = find(clInfo.ActiveUnit == true & clInfo.shank == 2);
s3 = find(clInfo.ActiveUnit == true & clInfo.shank == 3);
s4 = find(clInfo.ActiveUnit == true & clInfo.shank == 4);
s5 = find(clInfo.ActiveUnit == true & clInfo.shank == 5);
s6 = find(clInfo.ActiveUnit == true & clInfo.shank == 6);
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
S1 = [s1Control, s1CNO];
[p1, h1] = ranksum(s1Control, s1CNO);
S2 = [s2Control, s2CNO];
[p2, h2] = ranksum(s2Control, s2CNO);
S3 = [s3Control, s3CNO];
[p3, h3] = ranksum(s3Control, s3CNO);
S4 = [s4Control, s4CNO];
[p4, h4] = ranksum(s4Control, s4CNO);
S5 = [s5Control, s5CNO];
[p5, h5] = ranksum(s5Control, s5CNO);
S6 = [s6Control, s6CNO];
[p6, h6] = ranksum(s6Control, s6CNO);

figure; subplot(1, 6, 1);
boxplot(S1);

subplot(1, 6, 2);
boxplot(S2);

subplot(1, 6, 3);
boxplot(S3);

subplot(1, 6, 4);
boxplot(S4);

subplot(1, 6, 5);
boxplot(S5);

subplot(1, 6, 6);
boxplot(S6);

dS1 = ((median(s1CNO) - median(s1Control))/median(s1Control))*100
dS2 = ((median(s2CNO) - median(s2Control))/median(s2Control))*100
dS3 = ((median(s3CNO) - median(s3Control))/median(s3Control))*100
dS4 = ((median(s4CNO) - median(s4Control))/median(s4Control))*100
dS5 = ((median(s5CNO) - median(s5Control))/median(s5Control))*100
dS6 = ((median(s6CNO) - median(s6Control))/median(s6Control))*100