function PD_wav(subjects_mat, epoch_length)

% Script to calculate wavelet transform of epochs of length epoch_length
% centered around infusion time.

load(subjects_mat)

sampling_freq = 1000;

freqs = 1:200; no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 500]; no_bands = size(bands, 1);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    %% Calculating epochs, relative infusion time.
    
    basetime = basetimes(fo);
    
    base_index = basetime*sampling_freq;
    
    no_pre_epochs = floor(base_index/epoch_length);
    
    start_index = base_index - no_pre_epochs*epoch_length;
    
    no_post_epochs = floor(min(1500, size(PD_dec, 1)/sampling_freq - basetime)*sampling_freq/epoch_length);
    
    no_epochs = no_pre_epochs + no_post_epochs;
    
    clear Spec
    
    [epoch_no, t] = deal(nan(no_epochs, 1));
    
    for e = 1:no_epochs
        
        epoch_no(e) = e - no_pre_epochs - 1;
        
        epoch_start = start_index + (e - 1)*epoch_length + 1;
        
        epoch_end = start_index + e*epoch_length;
        
        t(e) = mean([epoch_start epoch_end])/sampling_freq - basetime;
        
        epoch_data = PD_dec(epoch_start:epoch_end, :);
        
        wavelet_spectrogram(epoch_data, sampling_freq, freqs, no_cycles, 1, [subj_name, '_', num2str(epoch_length/sampling_freq),'s_epoch_', num2str(e), 'wav.mat'])
        
    end
    
end