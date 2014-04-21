load('PD_dec.mat')

% sampling_rate=1000;

no_secs=floor(length(PD_dec)/sampling_rate);

t=[1:no_secs*sampling_rate]/sampling_rate;

% %% Morlet wavelet spectrogram.
% 
% MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
% freq = 1:200;
% scales = 1./(freq*MorletFourierFactor);
% 
% PD_sig =  struct('val',PD_dec,'period',1/sampling_rate);
% cwtPD = cwtft(PD_sig,'scales',scales);
% scales = cwtPD.scales;
% freq = 1./(scales*MorletFourierFactor);
% 
% figure;imagesc(t,freq,zscore(abs(cwtPD.cfs)));set(gca,'YDir','normal');
% xlabel('time (sec)'); ylabel('Pseudo-frequency');
% title('Morlet spectrogram of PD data')
% set(gca,'YTick',freq(1:10:length(freq)))
% grid on
% ylim([freq(1) freq(end)])

%% Gabor spectrogram.
L       = 10;
NFFT    = 2^L;
tw      = -NFFT/2+1:NFFT/2;
sigma   = .2;%[sec]
sigSamp = sigma*sampling_rate;
w       = sqrt(sqrt(2)/sigSamp)*exp(-pi*tw.*tw/sigSamp/sigSamp);
overlap = NFFT-1;
[SpecPD, T, F]=spectrogram(PD_dec,w,overlap,NFFT,sampling_rate);%deriv gaussian windowed spectrogram
figure;imagesc(F,T,zscore(abs(SpecPD'))');set(gca,'YDir','normal');
xlabel('time (sec)'); ylabel('Hz');
title('Gabor spectrogram of PD data')
set(gca,'YTick',freq(1:10:length(freq)))
grid on
ylim([freq(1) freq(end)])