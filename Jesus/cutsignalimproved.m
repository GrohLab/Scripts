[a1 a2 a3 a4]=pulsewidth(piezo);%obtain flags, test is a portion of a signal as an example.a2 contains the onset flags,
%a3 contains the offset flags. a4 contains the time duraction between onset
%and offset in time.
%fs = 1/signal sample interval
%make sure fs is equal in all signals!! fs=fs/sampling interval
a5=round(a2); 
fs=1/Control_LTP_Ch68.interval;
t=[1/fs:1/fs:size(v)/fs];%time scale 
%onset flag
%v1 is just a lfp signal to be plotted loaded in the workspace colum vector
ct=[];
for i=1:size(a5,1);
    spike=Control1(a5(i)-round(0.001*fs):1:a5(i)+round(0.1*fs));%cut signal
    ct(:,i)=spike;
end
a6=0;
figure
   subplot(3,1,1);
    plot(ct)
  wavecontrol=mean(ct(:,1:150),2);
  waveafterind=mean(ct(:,151:390),2);
  subplot(3,1,2)
 plot(wavecontrol);
 hold on;
 plot(waveafterind);
 hold off
 subplot(3,1,3)
 plot(mean(ct,2));
 allwave=mean(ct,2);
 hold off
 figure; plot(allwave);