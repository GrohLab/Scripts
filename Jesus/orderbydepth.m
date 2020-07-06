getConditionSignalsBF(fopen('Jittering_3mW_2.smrx'));
getDelayProtocol(pwd)
load('Jittering_3mW_2analysis.mat', 'Triggers')
Triggers
fs
lfp_spectrum = fftshift(fft(Triggers.lfp));
frequency_axis = -(fs/2):fs/length(lfp_spectrum):(fs/2) - fs/length(lfp_spectrum);
length(frequency_axis)
length(Triggers.lfp)
length(lfp_spectrum)
fs
fs/length(lfp_spectrum)
figure; plot(frequency_axis, lfp_spectrum)
figure; plot(frequency_axis, 10*log10(abs(lfp_spectrum)))
fs
s = spectrogram(Triggers.lfp, round(fs), round(fs * 0.1));
size(s)
spectrogram_timeaxis = (0:size(s,1)-1)/fs;
fs/2
fs/size(s,2)
fs/(size(s,2)*2)
figure; imagesc(20*log10(abs(s)))
spectrogram_frequencyaxis = (0:size(s,2)-1)*(fs/size(s,2));
max(spectrogram_frequencyaxis)
spectrogram_frequencyaxis = (0:size(s,2)-1)*(fs/(2*size(s,2)));
max(spectrogram_frequencyaxis)
timeaxis = (0:length(Triggers.lfp)-1)/fs;
max(timeaxis)
max(spectrogram_timeaxis)
spectrogram_timeaxis = (0:size(s,1)-1)/fs
clc
1/fs
size(s)
size(s)/fs
figure; plot(timeaxis, Triggers.lfp)
total_duration = max(timeaxis)
spectrogram_timeaxis = 0:size(s,1)/total_duration:total_duration-1;
max(spectrogram_timeaxis)
size(s,1)/total_duration
total_duration/(size(s,1)/total_duration)
spectrogram_timeaxis = 0:total_duration/size(s,1):total_duration-1;
max(spectrogram_timeaxis)
spectrogram_timeaxis = 0:total_duration/size(s,1):total_duration;
max(spectrogram_timeaxis)
spectrogram_timeaxis = 0:total_duration/(size(s,1)-1):total_duration;
max(total_duration)
total_duration/(size(s,1)-1)
size(spectrogram_timeaxis)
fs/2
1/(total_duration/(size(s,1)-1))
spectrogram_timeaxis
clc
size(s)
figure; imagesc(20*log10(abs(s)))
spectrogram_timeaxis = 0:total_duration/(size(s,2)-1):total_duration;
max(spectrogram_timeaxis)
close all
fs/size(s,1)
delta_frequency = fs/size(s,1);
spectrogram_frequencyaxis = 0:fs/(size(s,1)-1):size(s,1);
fs/2
max(s,1)
max(s)
spectrogram_frequencyaxis = 0:(fs/2)/size(s,1):fs/2;
size(spectrogram_frequencyaxis)
spectrogram_frequencyaxis = 0:(fs/2)/size(s,1):fs/2-1;
size(spectrogram_frequencyaxis)
spectrogram_frequencyaxis = 0:(fs/2)/size(s,1):(fs/2)-1;
size(spectrogram_frequencyaxis)
spectrogram_frequencyaxis = 0:(fs/2)/size(s,1):(fs/2);
size(spectrogram_frequencyaxis)
fs/2
fs/2-1
spectrogram_frequencyaxis = (0:(fs/2)/size(s,1):(fs/2))-1;
spectrogram_frequencyaxis = 0:(fs/2)/size(s,1)-1:(fs/2);
size(spectrogram_frequencyaxis)
spectrogram_frequencyaxis = 0:(fs/2)/(size(s,1)-1):(fs/2);
size(spectrogram_frequencyaxis)
figure; imagesc(spectrogram_timeaxis,spectrogram_frequencyaxis,20*log(abs(s))))
figure; imagesc(spectrogram_timeaxis,spectrogram_frequencyaxis,20*log(abs(s)))
figure; imagesc(spectrogram_timeaxis,spectrogram_frequencyaxis,100*log(abs(s)))
clear