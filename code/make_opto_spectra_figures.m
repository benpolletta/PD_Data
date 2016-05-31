function make_opto_spectra_figures(group_prefix, peak_suffix)

load([group_prefix, '_subjects.mat'], 'chan_labels', 'pd_labels')

no_chans = length(chan_labels);

no_pds_plotted = 2;

figure

for ch = 1:no_chans
    
    load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_0-200Hz_10trials_pct_spectrum_ch', num2str(ch), '_data_for_plot.mat'])

    All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se(:, 1:no_pds_plotted);
    
    [sig_lower, sig_higher] = find_sig(All_mean_mean(:, 1:no_pds_plotted), All_mean_ci);
    
    %% Plotting broadband spectrum.
    
    subplot(2, no_chans, ch)
    
    boundedline((1:200)', All_mean_mean(:, 1:no_pds_plotted), prep_for_boundedline(All_mean_ci), 'cmap', [0 1 1; 1 0 1]);
    
    axis tight
    
    add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
    
    add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
    
    set(gca, 'FontSize', 16)
    
    title(chan_labels{ch}, 'FontSize', 20)
    
    if ch == 1
        
        ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Pre-laser)'}, 'FontSize', 16)
        
    else
        
        legend(pd_labels(1:no_pds_plotted))
        
    end
    
    %% Plotting zoomed spectra.
    
    subplot(2, no_chans, 2 + ch)
    
    boundedline((1:50)', All_mean_mean(1:50, 1:no_pds_plotted), prep_for_boundedline(All_mean_ci(1:50, :)), 'cmap', [0 1 1; 1 0 1])
    
    axis tight
    
    add_stars(gca, (1:50)', logical(sig_lower(1:50)), 0, [1 .5 0])
    
    add_stars(gca, (1:50)', logical(sig_higher(1:50)), 1, [1 0 0])
    
    % plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
    
    set(gca, 'FontSize', 16)
    
    xlabel('Freq. (Hz)', 'FontSize', 16)
    
    if ch == 1
        
        ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Pre-laser)'}, 'FontSize', 16)
        
    end
    
end

% %% Plotting narrowand spectrum.
%     
% for ch = 1:no_chans
% 
%     load([group_prefix, '_8-30Hz_3-7cycles_6bands_8-30Hz_10trials_pct_spectrum_ch1_data_for_plot.mat'])
%     
%     subplot(3, no_chans, no_chans + ch)
%     
%     boundedline((8:.5:30)', WT_mean(:, 1:no_pds_plotted), prep_for_boundedline(WT_ci(:, 1:no_pds_plotted)))
%     
%     axis tight
%     
%     set(gca, 'FontSize', 16)
%     
%     % xlabel('Freq. (Hz)', 'FontSize', 16)
%     
%     if ch == 1
%         
%         ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)
%         
%     end
%     
% end

save_as_eps(gcf, [group_prefix, peak_suffix, '_spectra'])

save_as_pdf(gcf, [group_prefix, peak_suffix, '_spectra'])

end
