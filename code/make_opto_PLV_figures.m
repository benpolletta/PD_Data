function make_opto_PLV_figures(group_prefix, peak_suffix)

load([group_prefix, '_subjects.mat'], 'chan_labels', 'pd_labels')

no_chans = length(chan_labels);

no_pds_plotted = 2;

figure

%% Plotting zoomed magnitude of coherence.

load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_0-200Hz_10trials_PLV_Coh_sec_data_for_plot.mat'])

[sig_higher, sig_lower] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));

subplot(2, 1, 1)

boundedline((1:50)', PLV_mean(1:50, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:50, 1:no_pds_plotted)))

add_stars(gca, (1:50)', logical(sig_higher(1:50) + sig_lower(1:50)), 1, [1 0 0])

hold on

plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

axis tight

set(gca, 'FontSize', 16)

% title(chan_labels{1}, 'FontSize', 20)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Coherence'}, 'FontSize', 16)

%% Plotting zoomed angle of coherence.

load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_0-200Hz_10trials_PLV_dP_sec_data_for_plot.mat'])

[sig_higher, sig_lower] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));

subplot(2, 1, 2)

boundedline((1:50)', (180/pi)*PLV_mean(1:50, 1:no_pds_plotted), prep_for_boundedline((180/pi)*PLV_ci(1:50, 1:no_pds_plotted)))

add_stars(gca, (1:50)', logical(sig_higher(1:50) + sig_lower(1:50)), 1, [1 0 0])

hold on

plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

plot((1:50)', zeros(1, 50), ':k')

axis tight

set(gca, 'FontSize', 16)

% title(chan_labels{1}, 'FontSize', 20)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Angle of Coherence'}, 'FontSize', 16)

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

save_as_eps(gcf, [group_prefix, peak_suffix, '_PLV'])

save_as_pdf(gcf, [group_prefix, peak_suffix, '_PLV'])

end
