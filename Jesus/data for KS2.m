vectorBSandM1CVG=zeros(13,max(round(mean(BSandM1CVG,2))))
for i = 1:size(BSandM1CVG,1)
 approxMeanBSandM1=round(mean(BSandM1CVG,2))
 repBSandM1=repmat(i,approxMeanBSandM1(i),1)
    vectorBSandM1CVG(i,1:size(repBSandM1))=repBSandM1
end

vectorBSandS1CVG=zeros(13,max(round(mean(BSandS1CVG,2))))
for i = 1:size(BSandS1CVG,1)
 approxMeanBSandS1=round(mean(BSandS1CVG,2))
 repBSandS1=repmat(i,approxMeanBSandS1(i),1)
    vectorBSandS1CVG(i,1:size(repBSandS1))=repBSandS1
end


%V=[0 1 0 ;0 0 0; 1 1 1 ]
%test1=zeros(719,591);

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
test2=test1+test+test3+test4
imagesc(test2)
colormapJes=[1 1 1; 0 0 0; 0 0 0;0 0 0]
    