function make_spectra_figures(group_prefix, peak_suffix)

load([group_prefix, '_subjects.mat'], 'chan_labels')

figure

%% Plotting broadband spectrum.

load([group_prefix, peak_suffix, '_pct_8-30Hz_high_2.5_min_secs_pct_spectrum_ch1_data_for_plot.mat'])

subplot(2, 1, 1)

boundedline((1:200)', WT_mean, prep_for_boundedline(WT_ci))

axis tight

set(gca, 'FontSize', 16)

title(chan_labels{1}, 'FontSize', 20)

% xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

%% Plotting zoomed spectra.

subplot(2, 1, 2)

boundedline((1:50)', WT_mean(1:50, :), prep_for_boundedline(WT_ci(1:50, :)))

hold on

plot([15 30; 15 30], repmat([all_dimensions(@min, WT_mean(1:50, :) - WT_ci(1:50, :)); all_dimensions(@max, WT_mean(1:50, :) + WT_ci(1:50, :))], 1, 2), '--r')

axis tight

set(gca, 'FontSize', 16)

% title(chan_labels{1}, 'FontSize', 20)

xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

% %% Plotting narrowand spectra.
% 
% spec_labels = {'8-13Hz', '13-18Hz', '18-25Hz'};
% 
% for s = 1:length(spec_labels)
%     
%     load([group_prefix, '_8-30Hz_3-7cycles_6bands_pct_', spec_labels{s}, '_high_2.5_min_secs_pct_spectrum_ch1_data_for_plot.mat'])
%     
%     subplot(3, 3, 3 + s)
%     
%     boundedline((8:.5:30)', WT_mean, prep_for_boundedline(WT_ci))
%     
%     axis tight
%     
%     set(gca, 'FontSize', 16)
%     
%     xlabel('Freq. (Hz)', 'FontSize', 16)
%     
%     if s == 1
%         
%         ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)
%         
%     end
%     
% end

save_as_eps(gcf, [group_prefix, peak_suffix, '_spectra'])

save_as_pdf(gcf, [group_prefix, peak_suffix, '_spectra'])

end