function make_PLV_figures(group_prefix, peak_suffix, folder_flag)

load([group_prefix, '_subjects.mat'], 'chan_labels', 'pd_labels')

no_chans = length(chan_labels);

no_pds_plotted = 2;

figure

%% Plotting zoomed magnitude of coherence.

load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV', folder_flag, '_Coh_sec_pct_data_for_plot.mat'])

PLV_mean = All_mean_mean; PLV_ci = All_mean_ci;

[sig_lower, sig_higher] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));

subplot(2, 1, 1)

boundedline((1:50)', PLV_mean(1:50, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:50, 1:no_pds_plotted)))

add_stars(gca, (1:50)', logical(sig_lower(1:50)), 0, [1 .5 0])

add_stars(gca, (1:50)', logical(sig_higher(1:50)), 1, [1 0 0])

hold on

% plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

plot((1:50)', zeros(1, 50), ':k')

axis tight

set(gca, 'FontSize', 16)

% title(chan_labels{1}, 'FontSize', 20)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Phase-Locking Magnitude (%\Delta BL)'}, 'FontSize', 16)

%% Plotting zoomed angle of coherence.

load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_PLV', folder_flag, '_dP_sec_data_for_plot.mat'])

PLV_mean = All_mean_mean; PLV_ci = All_mean_ci;

[sig_lower, sig_higher] = find_sig(PLV_mean(:, 1:2), PLV_ci(:, 1:2));

subplot(2, 1, 2)

boundedline((1:50)', PLV_mean(1:50, 1:no_pds_plotted), prep_for_boundedline(PLV_ci(1:50, 1:no_pds_plotted)))

add_stars(gca, (1:50)', logical(sig_lower(1:50)), 0, [1 .5 0])

add_stars(gca, (1:50)', logical(sig_higher(1:50)), 1, [1 0 0])

hold on

% plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

plot((1:50)', zeros(1, 50), ':k')

axis tight

set(gca, 'FontSize', 16)

% title(chan_labels{1}, 'FontSize', 20)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Phase-Locking Angle'; '(Striatum - Motor Cortex)'}, 'FontSize', 16)

save_as_eps(gcf, [group_prefix, peak_suffix, '_PLV'])

save_as_pdf(gcf, [group_prefix, peak_suffix, '_PLV'])

end
