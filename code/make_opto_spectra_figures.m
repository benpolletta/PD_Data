function make_opto_spectra_figures(group_prefix, peak_suffix)

load([group_prefix, '_subjects.mat'], 'chan_labels', 'pd_labels')

no_chans = length(chan_labels);

no_pds_plotted = 2;

figure

for ch = 1:no_chans
    
    load([group_prefix, peak_suffix, '_pct_0-200Hz_10trials_pct_spectrum_ch', num2str(ch), '_data_for_plot.mat'])

    [sig_higher, sig_lower] = find_sig(WT_mean(:, 1:2), WT_ci(:, 1:2));
    
    %% Plotting broadband spectrum.
    
    subplot(2, no_chans, ch)
    
    boundedline((1:200)', WT_mean(:, 1:no_pds_plotted), prep_for_boundedline(WT_ci(:, 1:no_pds_plotted)))

    add_stars(gca, (1:200)', logical(sig_higher + sig_lower), 1, [1 0 0])
    
    axis tight
    
    set(gca, 'FontSize', 16)
    
    title(chan_labels{ch}, 'FontSize', 20)
    
    if ch == 1
        
        ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Pre-laser)'}, 'FontSize', 16)
        
    else
        
        legend(pd_labels(1:no_pds_plotted))
        
    end
    
    %% Plotting zoomed spectra.
    
    subplot(2, no_chans, 2 + ch)
    
    boundedline((1:50)', WT_mean(1:50, 1:no_pds_plotted), prep_for_boundedline(WT_ci(1:50, 1:no_pds_plotted)))

    add_stars(gca, (1:50)', logical(sig_higher(1:50) + sig_lower(1:50)), 1, [1 0 0])
    
    hold on
    
    plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')
    
    axis tight
    
    set(gca, 'FontSize', 16)
    
    % title(chan_labels{1}, 'FontSize', 20)
    
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

save_as_eps(gcf, [group_prefix, '_spectra'])

save_as_pdf(gcf, [group_prefix, '_spectra'])

end
