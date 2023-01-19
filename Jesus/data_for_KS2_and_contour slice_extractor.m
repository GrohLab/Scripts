vectorBSandM1CVG=zeros(13,max(round(mean(BSandM1CVG,2))))
for i = 1:size(BSandM1CVG,1)
 approxMeanBSandM1=round(mean(BSandM1CVG,2))
 repBSandM1=repmat(i,approxMeanBSandM1(i),1)
    vectorBSandM1CVG(i,1:size(repBSandM1))=repBSandM1
end
vectorBSandM1CVG(vectorBSandM1CVG==0)=nan

vectorBSandS1CVG=zeros(13,max(round(mean(BSandS1CVG,2))))
for i = 1:size(BSandS1CVG,1)
 approxMeanBSandS1=round(mean(BSandS1CVG,2))
 repBSandS1=repmat(i,approxMeanBSandS1(i),1)
    vectorBSandS1CVG(i,1:size(repBSandS1))=repBSandS1
end

vectorBSandS1CVG(vectorBSandS1CVG==0)=nan


vectorS1andM1CVG=zeros(13,max(round(mean(M1andS1CVG,2))))
for i = 1:size(M1andS1CVG,1)
 approxMeanS1andM1=round(mean(M1andS1CVG,2))
 repS1andM1=repmat(i,approxMeanS1andM1(i),1)
    vectorS1andM1CVG(i,1:size(repS1andM1))=repS1andM1
end
vectorS1andM1CVG(vectorS1andM1CVG==0)=nan

%%
vectorBSRC=zeros(13,max(round(mean(BSRC,2))))
for i = 1:size(BSRC,1)
 approxMeanBSRC=round(mean(BSRC,2))
 repBSRC=repmat(i,approxMeanBSRC(i),1)
    vectorBSRC(i,1:size(repBSRC))=repBSRC
end
vectorBSRC(vectorBSRC==0)=nan
vectorM1RC=zeros(13,max(round(mean(M1RC,2))))
for i = 1:size(M1RC,1)
 approxMeanM1RC=round(mean(M1RC,2))
 repM1RC=repmat(i,approxMeanM1RC(i),1)
    vectorM1RC(i,1:size(repM1RC))=repM1RC
end
vectorM1RC(vectorM1RC==0)=nan
vectorS1RC=zeros(13,max(round(mean(S1RC,2))))
for i = 1:size(S1RC,1)
 approxMeanS1RC=round(mean(S1RC,2))
 repS1RC=repmat(i,approxMeanS1RC(i),1)
    vectorS1RC(i,1:size(repS1RC))=repS1RC
end
vectorS1RC(vectorS1RC==0)=nan


i

vectorM1RCT=zeros(13,max(sum(M1RC,2)))
for i = 1:size(M1RC,1)
 approxMeanM1RC=sum(M1RC,2)
 repM1RC=repmat(i,approxMeanM1RC(i),1)
    vectorM1RCT(i,1:size(repM1RC))=repM1RC
end
vectorM1RCT(vectorM1RCT==0)=nan


%% extract contours slices
%V=[0 1 0 ;0 0 0; 1 1 1 ]
%close all
V=T(:,:,1)<=240;%grey values on inverted jpeg file
V=imresize(V,0.1);
test1=[]
for i=1:size(V,1);
    for ii=1:size(V,2);
        if V(i,ii)==0;
            test1(i,ii)=0;
        elseif V(i,ii)==1;
            test1(i,ii)=200;
            break
        end
    end
end

 
%V=[1 0 1 0 ;0 0 0 1; 1 1 1 1;0 1 0 1 ];
%test=zeros(4,4);
test=[]
for ii=1:size(V,2);
     for i=1:size(V,1);
        if V(i,ii)==0;
            test(i,ii)=0;
        elseif V(i,ii)==1;
            test(i,ii)=200;
            break
        end
    end
end

test3=[]

for ii=flip(1:size(V,2));
     for i=flip(1:size(V,1));
        if V(i,ii)==0;
            test3(i,ii)=0;
        elseif V(i,ii)==1;
            test3(i,ii)=200;
            break
        end
    end
end   

test4=[]

for i=flip(1:size(V,1));
    for ii=flip(1:size(V,2));
        if V(i,ii)==0;
            test4(i,ii)=0;
        elseif V(i,ii)==1;
            test4(i,ii)=200;
            break
        end
    end
end
test2=test1+test%+test3+test4
figure
imshow(T);
h=gca;h.Visible='on'
figure
imagesc(test2)
colormapJes=[1 1 1; 0 0 0; 0 0 0;0 0 0]
colormap(colormapJes); axis ij
[x y]=find(test2>=1);
figure
plot(y/61,x/61,'k*','MarkerSize',1); axis ij
hold on


%change values of V for labeled cell area detection 
V1=T(:,:,1)<=70;
V1=imresize(V1,0.1);
[x1,y1]=find(V1==1);
plot(y1/61,x1/61,'ro','MarkerSize',3,'MarkerFaceColor',[1 0 0]); axis ij %pixel size is 0.61*factor of 10 of IMresize
k=boundary(x1,y1,1);
patch(y1(k)/61,x1(k)/61,[1, 0 0],'FaceAlpha',0.5,'LineStyle','none')

%%
figure
plot(NeuronLoc(1:291,3),NeuronLoc(1:291,4),'ro','MarkerFaceColor','r')
hold on
plot(NeuronLoc(292:end,3),NeuronLoc(292:end,4),'go','MarkerFaceColor','g')
%%
figure

plot((xmc/0.61)*10,(ymc/0.61)*10,'r*','MarkerSize',0.7,'MarkerEdgeColor',[1 0 0])
hold on
kmc=boundary(x1mc,y1mc,0.9)
p1=patch((x1mc(kmc)/0.61)*10,(y1mc(kmc)/0.61)*10,[1 0 0],'FaceAlpha',0.4,'LineStyle','none')
plot((xsc/0.61)*10,(ysc/0.61)*10,'b*','MarkerSize',0.7,'MarkerEdgeColor',[0 0 1])
ksc=boundary(x1sc,y1sc,0.9)
p2=patch((x1sc(ksc)/0.61)*10,(y1sc(ksc)/0.61)*10,[0 0 1],'FaceAlpha',0.4,'LineStyle','none')
plot((xbs/0.61)*10,(ybs/0.61)*10,'g*','MarkerSize',0.7,'MarkerEdgeColor',[0 1 0])
kbs=boundary(x1bs,y1bs,0.9)
p3=patch((x1bs(kbs)/0.61)*10,(y1bs(kbs)/0.61)*10,[0 1 0],'FaceAlpha',0.4,'LineStyle','none')
title('Location of MC, BC, BS recepient areas in SC ')
legend([p1 p2 p3],{'MC','BC','BS'},'Location','southeast')
set(gca,'XTick',0:500:2500)
%set(gca,'YTickLabel',["0","0.5","1","1.5","2.0","2.5"
ylim([0 2800])
xlim([0 2800])
view([90 90]);
xlabel('Ventral<---->Dorsal (µm)')
ylabel('Lateral<---->Medial (µm)')