a=yControlOn-yL200On
b=[]
for i=1:size(a,2)
    if a(:,i)>0
        p=a(i);
    else b=a(i);
        ShorterMedianClusters(:,i)=b;
        LongerMedianClusters(:,i)=p;
        
    end
end