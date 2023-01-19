T=imfinfo('D:\C2-TestT.tif');
T1=T.ImageDescription
numbers = str2double(regexp(T1, '(?<=slices[^0-9])[0-9][0-9]','match'))
for i=1:numbers
T=imread('D:\C2-TestT.tif',i)*10;
testM(:,:,i)=T;
end
V=imresize(testM,0.25);
V1=subvolume(V,[150,270,250,350,10 20])