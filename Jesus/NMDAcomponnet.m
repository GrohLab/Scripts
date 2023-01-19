%% NMDA COMPONENT DISECTION
posMod=MIevok>=0.0001;negMod=MIevok<=-0.0001

PotS1C=nan(68,1000)
 spikesJES=(relativeSpkTmsStruct(1).SpikeTimes(wruIdx&posMod&signMod,1:end))
for i= 1:size(spikesJES,1)
   Test1=[(spikesJES{i,1:end})]
   PotS1C(i,1:size(Test1,2))=Test1
end
PotS1AI1=nan(68,1000)
 spikesJES1=(relativeSpkTmsStruct(2).SpikeTimes(wruIdx&posMod&signMod,1:150))
for i= 1:size(spikesJES1,1)
   Test2=[(spikesJES1{i,:})]
   PotS1AI1(i,1:size(Test2,2))=Test2
end

PotS1AI2=nan(68,1000)
 spikesJES2=(relativeSpkTmsStruct(2).SpikeTimes(wruIdx&posMod&signMod,150:300))
for i= 1:size(spikesJES2,1)
   Test3=[(spikesJES2{i,1:end})]
   PotS1AI2(i,1:size(Test3,2))=Test3
end

kernMatrix =nan(68,1000);figure
subplot(1,2,1)
for i=1:size(PotS1C,1)
[yCtrE xCtrE]=ksdensity(PotS1C(i,:),'Bandwidth',0.0005,'Kernel','epanechnikov','NumPoints',1000)
kernMatrix(i,1:1000)=yCtrE;
plot(xCtrE,kernMatrix(i,:)); hold on
end
subplot(1,2,2)
imagesc([-37:0.12:83-0.12],[1:1:68],kernMatrix)
colormap(jet(128))
camlight; lighting gouraud

kernMatrixAI =nan(68,240);figure
subplot(1,2,1)
for i=1:size(PotS1AI1,1)
%[yCtrE xCtrE]=ksdensity(PotS1AI1(i,:),'Bandwidth',0.0005,'Kernel','epanechnikov','NumPoints',1000)
%kernMatrixAI(i,1:1000)=yCtrE;
yAI1=histcounts(PotS1AI1(i,:),240);
plot([-37:0.5:83-0.12],yAI1); hold on
kernMatrixAI(i,1:240)=yAI1
end
subplot(1,2,2)
imagesc([-37:0.12:83-0.12],[1:1:68],kernMatrixAI)
colormap(jet(128))


kernMatrixAI2 =nan(68,1000);figure
subplot(1,2,1)
for i=1:size(PotS1AI2,1)
[yCtrE xCtrE]=ksdensity(PotS1AI2(i,:),'Bandwidth',0.0005,'Kernel','epanechnikov','NumPoints',1000)
kernMatrixAI2(i,1:1000)=yCtrE;
plot(xCtrE,kernMatrixAI2(i,:)); hold on
end
subplot(1,2,2)
imagesc([-37:0.12:83-0.12],[1:1:68],kernMatrixAI2)
colormap(jet(128))