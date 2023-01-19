
[num,txt,raw] = xlsread('Z:\Jesus\TestCellExplorerSCBehavior\211205_C\Behaviour\roller2021-12-05T14_28_44DLC_resnet50_AwakenSCJul20shuffle1_1030000.xlsx','K:M')%ReadExcelData
puffTimes=Conditions(1).Triggers %puff in samples
PuffTimes1=puffTimes./fs % puff Times in time using fs of the intan recording 30000hz
puffTimes2=puffTimes*cfs %stimation of the position sample in the camera recording fs
zsczz=[];
for i=1:size(spikes.times,2) % all clusters aligned with camera fs
cfs=653.98
Fr1=spikes.times{1,i};%firing rate primera celula optotag
Fr2=diff(Fr1);
Fr3=1./Fr2;
Fr3=[Fr3;0];
AsocPos1=[];
AsocPos1=zeros(size(num(:,1)));%preallocate
equivPos1=Fr1*cfs; % equivalent position (Time*fs camera
equivPos1=round(equivPos1);
AsocPos1(equivPos1)=Fr3;
xx=[deeplabcut.position.x deeplabcut.position.x];
yy=[deeplabcut.position.y deeplabcut.position.y];
zz=[AsocPos1 AsocPos1];
zsczz(1:size(AsocPos1,1),i)=AsocPos1;%zscore(AsocPos1);
AsocPosAll(:,i)=AsocPos1;
sp=[deeplabcut.speed deeplabcut.speed];
%idj1=AsocPos1>2;
end
figure
plot(deeplabcut.timestamps,deeplabcut.position.y)
hold on
yyaxis right
%plot(deeplabcut.timestamps(idj1),zz(idj1),'*r')
plot(deeplabcut.timestamps,zscore(zz),'*r')

%[a1 a2 a3 a4]=pulsewidth(test);%obtain flags, test is a portion of a signal as an example.a2 contains the onset flags,
%a3 contains the offset flags. a4 contains the time duraction between onset
%and offset in time.
%fs = 1/signal sample interval
a5=round(a2); 
t=[1/fs:1/fs:size(v)/fs];%time scale 
%onset flag
%v1 is just a lfp signal to be plotted loaded in the workspace colum vector
ct=[];
for i=1:size(puffTimes2',2); %roller speed
    spike=vels(puffTimes2(i)-500:1:puffTimes2(i)+500);%cut signal
    ct(:,i)=spike;
end

ctAsocPosAll=[]; %firing rate neurons 3D matrix timexUIDxtrial
for i=1:size(puffTimes2',2);
    for ii=1:size(AsocPosAll,2);
    AsocPosAlli=AsocPosAll(:,ii);
    spike1=AsocPosAlli(puffTimes2(i)-500:1:puffTimes2(i)+500);%cut signal
    ctAsocPosAll(:,ii,i)=spike1;
    sprintf("cycle1 %d\n",ii)
    end
    sprintf("cycle2 %d\n",i)
end
cutAllclusPuff=mean(ctAsocPosAll,3); %mean per trial

delPufAsocPosAll=AsocPosAll;
for i=1:size(puffTimes2',2)
    tobedeletT=(puffTimes2(i)-500:1:puffTimes2(i)+500)
    tobedeleted(:,i)=tobedeletT;
end
iDdel=round(tobedeleted(:)) %Id positions of the array to be deleted
delPufAsocPosAll(iDdel,:)=[]; 

a=deeplabcut.position.x;%whisking in x
awp=a; 
awp(iDdel)=[];%whisking in x without -0.5 to 0.5 airpuff period.
b=deeplabcut.position.y; %whisking in y
bwp=b;
bwp(iDdel)=[]; %whisking in y without -0.5 to 0.5 airpuff period.

%Evoked whsiking in y 
EvokedB=[];
for i=1:size(puffTimes2',2); %roller speed vs y-whisking
    EvokedT=b(puffTimes2(i)-500:1:puffTimes2(i)+500);%cut signal
    EvokedB(:,i)=EvokedT;
end

EvokedA=[];
for i=1:size(puffTimes2',2); %roller speed vs x-whiksing
    EvokedTA=a(puffTimes2(i)-500:1:puffTimes2(i)+500);%cut signal
    EvokedA(:,i)=EvokedTA;
end

evokedN=[] %nose 
for i=1:size(puffTimes2',2); %roller speed
    EvokedTN=rw(puffTimes2(i)-500:1:puffTimes2(i)+500);%cut signal
    EvokedN(:,i)=EvokedTN;
end

%Laser Control without puff
LaserTimesC=Conditions(4).Triggers;
LaserTimes1C=LaserTimes./fs;
LaserTimes2C=LaserTimes1C.*cfs

ctAsocPosAllLas=[]; %firing rate neurons 3D matrix timexUIDxtrial
for i=1:size(LaserTimes2C',2);
    for ii=1:size(AsocPosAll,2);
    AsocPosAllLi=AsocPosAll(:,ii);
    spike2=AsocPosAllLi(LaserTimes2C(i)-500:1:LaserTimes2C(i)+500);%cut signal
    ctAsocPosAllLas(:,ii,i)=spike2;
    sprintf("cycle1 %d\n",ii)
    end
    sprintf("cycle2 %d\n",i)
end
cutAllclusLas=mean(ctAsocPosAllLas,3); %mean per trial
%only laser
EvokedBLC=[];
for i=1:size(LaserTimes2C',2); %only laser pulses
    EvokedTLC=b(LaserTimes2C(i)-500:1:LaserTimes2C(i)+500);%cut signal
    EvokedBLC(:,i)=EvokedTLC;
end



EvokedALC=[];
for i=1:size(LaserTimes2C',2); %only laser pulses
    EvokedTALC=a(LaserTimes2C(i)-500:1:LaserTimes2C(i)+500);%cut signal
    EvokedALC(:,i)=EvokedTALC;
end

%only puff

PuffTimesonlyC=Conditions(5).Triggers;
PuffTimesonly1C=PuffTimesonlyC./fs;
PuffTimesonly2C=PuffTimesonly1C.*cfs

ctAsocPosAllPonly=[]; %firing rate neurons 3D matrix timexUIDxtrial 
for i=1:size(PuffTimesonly2C',2);
    for ii=1:size(PuffTimesonly2C,2);
    AsocPosAllPi=AsocPosAll(:,ii);
    spike3=AsocPosAllPi(PuffTimesonly2C(i)-500:1:PuffTimesonly2C(i)+500);%cut signal
    ctAsocPosAllPonly(:,ii,i)=spike3;
    sprintf("cycle1 %d\n",ii)
    end
    sprintf("cycle2 %d\n",i)
end
cutAllclusPonly=mean(ctAsocPosAllPonly,3); %mean per trial
EvokedBPC=[];
for i=1:size(LaserTimes2C',2); %only puff
    EvokedTPC=b(LaserTimes2C(i)-500:1:LaserTimes2C(i)+500);%cut signal
    EvokedBPC(:,i)=EvokedTPC;
end

EvokedALC=[];
for i=1:size(LaserTimes2C',2); %only puff
    EvokedTAPC=a(LaserTimes2C(i)-500:1:LaserTimes2C(i)+500);%cut signal
    EvokedAPC(:,i)=EvokedTAPC;
end


%%
figure
hs=surf(xx,yy,zz,'EdgeColor','interp','FaceColor','interp')
caxis([0 200])
colormap(jet(64))
view(3)
caxis([0 100])
figure
hs=surf(yy,yy,zz,'EdgeColor','interp','FaceColor','interp')
caxis([0 200])
colormap(jet(64))
view(2)
caxis([0 100])
figure
hs=surf(xx,xx,zz,'EdgeColor','interp','FaceColor','interp')
caxis([0 200])
colormap(jet(64))
view(3)
caxis([0 100])