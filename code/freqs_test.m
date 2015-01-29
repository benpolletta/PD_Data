function [max_freqs, max_std_freqs] = freqs_test(folder, lim_mins, basetime, chan, freqs, cycles)

load(sprintf('%s/%s_all_channel_data_dec.mat', folder, folder([1:2,4:6])))

lims = (lim_mins*60 + basetime)*500;

PD_chosen = PD_dec(lims(1):lims(2), chan);

Spec_chosen = wavelet_spectrogram(PD_chosen, 500, freqs, cycles, 0, '');

[~, max_freq_indices] = max(Spec_chosen, [], 2);

max_freqs = freqs(max_freq_indices);

[~, max_std_freq_indices] = max(zscore(Spec_chosen), [], 2);

max_std_freqs = freqs(max_std_freq_indices);

for i = 1:((diff(lims, [], 2) + 1)/(10*500))
    
    indices = ((i - 1)*5000 + 1):(i*5000);

    t = (lims(1) + indices - basetime)/500;
    
    figure
    
    subplot(2,1,1)
    
    imagesc(t, freqs, abs(Spec_chosen(indices, :))')
    
    axis xy
    
    hold on
    
    plot(t, max_freqs(indices)', 'w--', 'LineWidth', 2)
    
    plot(t, max_std_freqs(indices)', 'y:', 'LineWidth', 2)
    
    axis xy, 
    
    subplot(2,1,2), 
    
    plot(t, PD_chosen(indices))
    
    axis tight
    
    box off

end