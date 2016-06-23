function make_spectra_figures(group_prefix, peak_suffix)

load([group_prefix, '_subjects.mat'], 'chan_labels')

figure

%% Plotting broadband spectrum.

load([group_prefix, '_1-200Hz_3-21cycles_7bands', peak_suffix, '_pct_15-30Hz_high_2.5_min_secs_pct_spectrum_missing_2_ch1_data_for_plot.mat'])

All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;

[sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);

subplot(2, 1, 1)

boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))

axis tight

add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])

add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])

set(gca, 'FontSize', 16)

title(chan_labels{1}, 'FontSize', 20)

ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

%% Plotting zoomed spectra.

subplot(2, 1, 2)

boundedline((1:50)', All_mean_mean(1:50, :), prep_for_boundedline(All_mean_ci(1:50, :)))

axis tight

add_stars(gca, (1:50)', logical(sig_lower(1:50)), 0, [1 .5 0])

add_stars(gca, (1:50)', logical(sig_higher(1:50)), 1, [1 0 0])

% plot([15 30; 15 30], repmat(ylim', 1, 2), '--r')

set(gca, 'FontSize', 16)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

save_as_eps(gcf, [group_prefix, peak_suffix, '_spectra'])

save_as_pdf(gcf, [group_prefix, peak_suffix, '_spectra'])

end