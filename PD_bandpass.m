load('PD_dec.mat')

no_secs=floor(length(PD_dec)/sampling_freq);

t=[1:no_secs*sampling_freq]/sampling_freq;

bands = [(110:180)'-10 (110:180)' (110:180)'+10];
bands = max(bands,0);

% [~,~,bands,S,H,A,P,~,~] = filter_fft_filtfilt(PD_dec,'sampling_freq',sampling_freq,'bands',bands);
[~,~,bands,S,H,A,P,~,~] = filter_fft_new(PD_dec,'sampling_freq',sampling_freq,'bands',bands);

figure()

imagesc(t,bands,zscore(A))
axis xy

save_as_pdf(gcf,'PD_spec_BW')
