[a1 a2 a3 a4]=pulsewidth(test);%obtain flags, test is a portion of a signal as an example.a2 contains the onset flags,
%a3 contains the offset flags. a4 contains the time duraction between onset
%and offset in time.
%fs = 1/signal sample interval
a5=round(a2); 
t=[1/fs:1/fs:size(v)/fs];%time scale 
%onset flag
%v1 is just a lfp signal to be plotted loaded in the workspace colum vector
ct=[];
for i=1:size(a5,1);
    spike=v1(a5(i)-50:1:a5(i)+30000);%cut signal
    ct(:,i)=spike
end
a6=0;
x2=[(a6-50:1:a6+30000)/fs];%time scale for evoked potentials
   subplot(2,1,1);
    plot(x2,ct,x2,test(a5));
  wave=mean(ct,2);
  subplot(2,1,2)
  plot(x2,wave)