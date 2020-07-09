%Control1 es una matriz que importo desde spike 2 a traves del script
%"raster dumper"

Jcontrol=control(:)
Jcontrol=Jcontrol*1000
ExcludeJcontrolOn=[find(Jcontrol<2);find(Jcontrol>50)]
Jcontrol(ExcludeJcontrolOn)=nan
JuxtaControlOn=Jcontrol(~isnan(Jcontrol))
figure;hist(JuxtaControlOn,30)
xlim([0 30])
title('Juxtacelular recordings')
hold on
%Laser

JL300=L300(:)
JL300=JL300*1000
ExcludeJL300On=[find(JL300<2);find(JL300>30)]
JL300(ExcludeJL300On)=nan
JuxtaL300On=JL300(~isnan(JL300))
figure;hist(JuxtaL300On,30)
xlim([0 30])
title('Juxtacelular recordings')

%cumulative fraction
cdfplot(JuxtaControlOn)
hold on
cdfplot(JuxtaL300On)
xlabel('time_{1ms bin}')
ylabel ('Cumulative Fraction')
legend ("Control","Laser")

figure;
plot(cumsum(hist(JuxtaControlOn)))
hold on; 
plot(cumsum(hist(JuxtaL300On)))


OnJuxtalabels={'control','laser'}
juxtaALLon=nan(size(JuxtaL300On,1),2)
juxtaALLon(1:size(JuxtaControlOn,1),1)=[JuxtaControlOn]
juxtaALLon(1:size(JuxtaL300On,1),2)=[JuxtaL300On]
boxplot(juxtaALLon,OnJuxtalabels)

scatterhist(JuxtaControlOn,JuxtaL300On(1:size(JuxtaControlOn)),'NBins',60,'Color','k','Marker','*','Kernel','overlay','Direction','out')
hold on
 plotAdded(fitlm(JuxtaControlOn,JuxtaL300On(1:size(JuxtaControlOn))))
 hold on
 ContourJuxta=hist3(juxtaALLon,[60 60])
figure, contour(ContourJuxta')
 
