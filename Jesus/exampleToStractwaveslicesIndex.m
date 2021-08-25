ct=[]
for i=1:size(b,2);
    spike=a(b(i)-10:1:b(i)+5);
    ct(:,i)=spike
end
    
wave=mean(ct,2);
