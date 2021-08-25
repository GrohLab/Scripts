
clear all
%open('D:\mat to brainstorm\TestAllChannels1.mat')
path='C:\Users\jesus\Desktop\Test'
myfile='Test18.mat'
cd(path)
names=fieldnames(open(myfile));
goods=[];
for i=1:numel(names)
    if ~isempty(strfind(names{i},'Ch'))
        if isempty(strfind(names{i},'Ch31'))
            goods=[goods;i];
        end
    end
end
names=names(goods)

%get number sample points
firstchan=load(myfile,names{1})
eval(['nsamples=firstchan.' names{1} '.length;'])
nchan=numel(names);
%%
V=zeros(nsamples,nchan);  %%preallocate space

%get every individual channel
for i=1:nchan
    thischan=names{i}; %name of channel as string
    v=load(myfile,thischan); %load channel from file
    eval([' V(:,i)=v.' thischan '.values;']) %put channel values in to preallocated matrix
    clear v
   
end

outputfile=[myfile(1:end-4) '_matrix']
save(outputfile,'V')
save(outputfile,'V','-v7.3')