function [max_freqs, max_std_freqs, max_pct_freqs] = freqs_test_2(folder, lim_mins, basetime, chan)

beta_freqs = 8:30;

load(sprintf('%s/%s_all_channel_data_dec.mat', folder, folder([1:2, 4:6])))

load(sprintf('%s/%s_wt.mat', folder, folder([1:2, 4:6])))

load(sprintf('%s/%s_wt_pct.mat', folder, folder([1:2, 4:6])))

lims = (lim_mins*60 + basetime)*500;

PD_chosen = PD_dec(lims(1):lims(2), chan);

Spec_chosen = Spec(lims(1):lims(2), beta_freqs, chan); % wavelet_spectrogram(PD_chosen, 500, freqs, cycles, 0, '')

Spec_pct_chosen = Spec_pct(lims(1):lims(2), beta_freqs, chan);

[~, max_freq_indices] = max(Spec_chosen, [], 2);

max_freqs = beta_freqs(max_freq_indices);

[~, max_std_freq_indices] = max(zscore(Spec_chosen), [], 2);

max_std_freqs = beta_freqs(max_std_freq_indices);

[~, max_pct_freq_indices] = max(Spec_pct_chosen, [], 2);

max_pct_freqs = beta_freqs(max_pct_freq_indices);

no_figs = (diff(lims, [], 2) + 1)/(10*500);

for i = 1:no_figs
    
    indices = ((i - 1)*5000 + 1):(i*5000);

    t = (lims(1) + indices - basetime)/500;
    
    figure
    
    subplot(2,1,1)
    
    imagesc(t, beta_freqs, abs(Spec_chosen(indices, :))')
    
    axis xy
    
    hold on
    
    h(1) = plot(t, max_freqs(indices)', 'w'); % , 'LineWidth', 2)
    
    h(2) = plot(t, max_std_freqs(indices)', 'w:', 'LineWidth', 2);
    
    h(3) = plot(t, max_pct_freqs(indices)', 'Color', 'w', 'LineStyle', '--', 'LineWidth', 2);
    
    if i == 1 
        
        ledge_h = legend(h, {'Raw','z-Scored','% Baseline Change'});
        
        set(ledge_h, 'Color', 'b', 'TextColor', 'w')
    
    end
    
    axis xy, 
    
    subplot(2,1,2), 
    
    plot(t, PD_chosen(indices))
    
    axis tight
    
    box off
    
    save_as_pdf(gcf, sprintf('freqs_test_%s_min%d-%d_ch%d_fig%d', folder, lim_mins, chan, i))

end