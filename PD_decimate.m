function PD_decimate(prefix)

load([prefix,'_all_data.mat'])

% PD_data = [fliplr(PD_data(1:10*sampling_rate)) PD_data fliplr(PD_data(1:10*sampling_rate))];
% 
% PD_data = eegfilt(PD_data, sampling_rate, 0, 300);
% 
% PD_data = PD_data((5*sampling_rate + 1):(end - 5*sampling_rate));
% 
% PD_data = detrend(PD_data, 'linear');

PD_dec = decimate(PD_data,10);

PD_dec = decimate(PD_dec,2);

sampling_freq = sampling_freq/20;

save([prefix,'_all_data_dec.mat'],'PD_dec','sampling_freq')