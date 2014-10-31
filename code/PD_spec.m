load('PD_dec.mat')

% PD_dec=decimate(PD_data,10);
% PD_dec=decimate(PD_dec,2);

% PD_dec=detrend(PD_dec,'linear');

% sampling_freq=1000;

no_secs=floor(length(PD_dec)/sampling_freq);

t=(1:no_secs*sampling_freq)/sampling_freq;

% [first_pmtm,f]=pmtm(PD_dec(1:sampling_freq),[],[],sampling_freq);
% 
% f(f>200)=[];
% 
% pmtm_length=length(f);
% 
% spectrum=zeros(no_secs,pmtm_length);
% 
% spectrum(1,:)=first_pmtm(1:pmtm_length);
% 
% for s=1:no_secs
%     
%     spec_temp=pmtm(detrend(PD_dec([1:sampling_freq]+(s-1)*sampling_freq),'linear'),[],[],sampling_freq);
%     
%     spectrum(s,:)=spec_temp(1:pmtm_length);
%     
% end
% 
% % spec_mean=ones(size(spectrum))*diag(nanmean(spectrum));
% % spec_std=ones(size(spectrum))*diag(nanstd(spectrum));
% % spec_norm=(spectrum-spec_mean)./spec_std;
% 
% figure()
% 
% imagesc(t,f,zscore(spectrum))
% axis xy
% 
% save_as_pdf(gcf,'PD_spec_by_10s')

bands = [(1:.5:50)'-3 (1:.5:50)' (1:.5:50)'+3];
bands = max(bands,0);

[~,~,bands,S,H,A,P,~,~] = filter_fft_filtfilt(PD_dec,'sampling_freq',sampling_freq,'bands',bands);

figure()

imagesc(t,bands,zscore(A))
axis xy

save_as_pdf(gcf,'PD_spec_BW')
