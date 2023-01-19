%% extract contours slices
%V=[0 1 0 ;0 0 0; 1 1 1 ]
close all
myfiles=dir('Z:\Berin\3D reconstruction inputs Cx_BS\*green.jpg');

for j=1:numel(myfiles);
    
%T=imread('Z:\Berin\3D reconstruction inputs Cx_BS\Rostral0 clear outside green.jpg');
 T=imread(fullfile(['Z:\Berin\3D reconstruction inputs Cx_BS\',[myfiles(j).name]]));
T=rgb2gray(T);
T=255-T;
T=flip(T);
V=T<=252;%grey values on inverted jpeg file
%V=imresize(V,0.1);

test1=[];
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
test=[];
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

test3=[];

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

test4=[];

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
test2=test4;%+test3;%+test+test1;
test2(1,:)=0;
%figure
%imshow(T);
h=gca;h.Visible='on'
%figure
%imagesc(test2)
colormapJes=[1 1 1; 0 0 0; 0 0 0;0 0 0];
colormap(colormapJes); axis ij
[x y]=find(test2>=1);
figure
plot(y/61,x/61,'k*','MarkerSize',1); axis ij
hold on


%change values of V for labeled cell area detection 
V1=T<=225;
%V1=imresize(V1,0.1);
[x1,y1]=find(V1==1);
plot(y1/61,x1/61,'go','MarkerSize',3,'MarkerFaceColor',[0 1 0]); axis ij %pixel size is 0.61*factor of 10 of IMresize
k=boundary(x1,y1,0.8);
patch(y1(k)/61,x1(k)/61,[0 1 0],'FaceAlpha',0.5,'LineStyle','none')
hold on

myfilesred=dir('Z:\Berin\3D reconstruction inputs Cx_BS\*red.jpg');

%T=imread('Z:\Berin\3D reconstruction inputs Cx_BS\Rostral0 clear outside red.jpg');
T=imread(fullfile(['Z:\Berin\3D reconstruction inputs Cx_BS\',[myfilesred(j).name]]));
T=rgb2gray(T);
T=255-T;
T=flip(T);
V1=T<=205;
%V1=imresize(V1,0.1);
[x1,y1]=find(V1==1);
plot(y1/61,x1/61,'ro','MarkerSize',3,'MarkerFaceColor',[1 0 0]); axis ij %pixel size is 0.61*factor of 10 of IMresize
k=boundary(x1,y1,0.8);
patch(y1(k)/61,x1(k)/61,[1, 0 0],'FaceAlpha',0.5,'LineStyle','none');
box off
hold off
%figure
%imshow(T)
%x3d1=x; y3d1=y;
%x13d1=x1; y13d1=y1;
x3d=[];y3d=[];
x3d(1:size(x,1),j)=x
y3d(1:size(y,1),j)=y
x3d1(1:size(x1,1),j)=x1
y3d1(1:size(y1,1),j)=y1
end

plot3((zeros(size(x3d(:,1),1),1)+0),flip(y3d(:,1)),flip(x3d(:,1)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
hold on
plot3((zeros(size(x3d(:,2),1),1)+50),flip(y3d(:,2)),flip(x3d(:,2)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,3),1),1)+100),flip(y3d(:,3)),flip(x3d(:,3)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,4),1),1)+150),flip(y3d(:,4)),flip(x3d(:,4)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,5),1),1)+200),flip(y3d(:,5)),flip(x3d(:,5)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,6),1),1)+250),flip(y3d(:,6)),flip(x3d(:,6)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,7),1),1)+300),flip(y3d(:,7)),flip(x3d(:,7)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,8),1),1)+350),flip(y3d(:,8)),flip(x3d(:,8)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
plot3((zeros(size(x3d(:,9),1),1)+400),flip(y3d(:,9)),flip(x3d(:,9)),'LineStyle','none','MarkerSize',0.1,'Marker','diamond','MarkerFaceColor',[0.1 0.1 0.1])
box off
axis ij
xlabel("x")
ylabel("y")


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