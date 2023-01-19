T=imfinfo('C:\Users\jesus\Downloads\Brain_Tetrodes_overlay.tif');
T1=T.ImageDescription
numbers = str2double(regexp(T1, '(?<=slices[^0-9])[0-9][0-9]','match'))
testM=[];
for i=1:numbers
T=imread('C:\Users\jesus\Downloads\Brain_Tetrodes_overlay.tif',i);
testM(:,:,i)=T(:,:,1);
end

V=imresize(testM,0.125);
V1=subvolume(V,[150,270,250,350,10 20])
slice(V,[],[],20)
hold on
slice(V,[],[],[10])
slice(V,[],[],1)
zlim([-5 25])
shading flat
camlight;
lighting gouraud
isosurface(V,200)
isocaps(V,20)
camlight; lighting  gouraud
%%
for i= 1:numbers
    slice(V,[],[],i)
    hold on;
end
shading flat
camlight; lighting gouraud;