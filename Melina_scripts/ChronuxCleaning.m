


 
clear all
cd 'Z:\Melina\MC-1-2019\2019-10-24-freely-optoExp-0-and-70percent\MC-1-2019_20191024_optoExp-0%,70%_Kilosort'
load dataCell_ds



%%
%this should be eventually a function, once we figure out what we want
for i=1 %:numel(dataCell_ds)
    lfp=dataCell_ds{i}';
    ppms=1;  %this is assuming that the downsampling was always to 1 ms!!!! hard-coded in laser plotting at the moment
    dsf=30  %this too
    T=numel(lfp)/(1000*ppms); %how long is recording
    params.Fs=(1000*ppms)/dsf; % sampling frequency
    params.fpass=[0 params.Fs/2]; % band of frequencies to be kept
    params.fpass=[0 500]; % band of frequencies to be kept
    Ktapers=9;
    NW=(Ktapers+1)/2;
    params.tapers=[NW Ktapers]; % taper parameters
    params.pad=2; % pad factor for fft
    params.err=[2 0.05];
    params.trialave=1;
    movingwin=[1 0.05];
    segave=1;
    wintrig=[5*movingwin(1) 5*movingwin(1)];
    winseg=2*movingwin(1);
    
    % get rid of slow drift and line noise (50 Hz)
    dataDetrend=locdetrend(lfp,params.Fs,[5 0.05]);
    %[f] = delineTest(dataDetrend,params.Fs)
   % dataDetrend = delineSignal(dataDetrend,params.Fs,(1:15)*f, 4);  %what is 1:9??? and 4???
   % dataDetrend=rmlinesc(dataDetrend,params);
    
    %%
    % compute segmented spectrum
    NT=5*round(winseg*params.Fs);
    data1=dataDetrend(1:NT)';
    clean.LFP_ds=dataDetrend';
    figure
    [clean.S,clean.f,clean.varS,clean.C,clean.Serr]=mtspectrumsegc(data1,winseg,params,segave);
    plot(clean.f,10*log10(clean.S),clean.f,10*log10(clean.Serr(1,:)),clean.f,10*log10(clean.Serr(2,:)));
    clean.ds=dsf;
    clean.params=params;
    %subplot(212)
    %imagesc(f,f,C); axis xy;colorbar;
    
    
   % rmfield(clean,'LFP')
   % rmfield(clean,'filteredResponse')
   % POm_clean{i}=clean;
   
   
   
   
end


%%


  thisCell=POm{i};
    clean=POm{i};
    lfp=thisCell.LFP';
    T=numel(lfp)/(1000*thisCell.ppms)
    
    dsf=20;
    lfp=downsample(lfp,dsf);
    
    %chronux parameters
    params.Fs=(1000*thisCell.ppms)/dsf; % sampling frequency
    params.fpass=[0 params.Fs/2]; % band of frequencies to be kept
    params.fpass=[0 500]; % band of frequencies to be kept
    Ktapers=9;
    NW=(Ktapers+1)/2;
    params.tapers=[NW Ktapers]; % taper parameters
    
    params.pad=2; % pad factor for fft
    params.err=[2 0.05];
    params.trialave=1;
    movingwin=[1 0.05];
    segave=1;
    wintrig=[5*movingwin(1) 5*movingwin(1)];
    winseg=2*movingwin(1);
    
    % get rid of slow drift and line noise (50 Hz)
    dataDetrend=locdetrend(lfp,params.Fs,[.5 0.05]);
    [f] = delineTest(dataDetrend,params.Fs)
    delinedEEG = delineSignal(dataDetrend,params.Fs,(1:15)*f, 4);  %what is 1:9??? and 4???
    delinedEEG=rmlinesc(delinedEEG,params);
    
    % compute segmented spectrum
    NT=5*round(winseg*params.Fs);
    data1=delinedEEG(1:NT)';
    clean.LFP_ds=delinedEEG';
    figure
    [clean.S,clean.f,clean.varS,clean.C,clean.Serr]=mtspectrumsegc(data1,winseg,params,segave);
    plot(clean.f,10*log10(clean.S),clean.f,10*log10(clean.Serr(1,:)),clean.f,10*log10(clean.Serr(2,:)));
    clean.ds=dsf;
    clean.params=params;
    %subplot(212)
    %imagesc(f,f,C); axis xy;colorbar;
    
    
    rmfield(clean,'LFP')
    rmfield(clean,'filteredResponse')
    POm_clean{i}=clean;