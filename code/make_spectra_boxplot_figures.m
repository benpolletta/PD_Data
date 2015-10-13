function make_spectra_boxplot_figures(group_prefix)

load([group_prefix, '_subjects.mat'], 'chan_labels')

figure

%% Plotting broadband spectrum.

load([group_prefix, '_pct_8-30Hz_high_2.5_min_secs_pct_spectrum_ch1_data_for_plot.mat'])

subplot(3, 1, 1)

boundedline((1:200)', WT_mean, prep_for_boundedline(WT_ci))

axis tight

set(gca, 'FontSize', 16)

title(chan_labels{1}, 'FontSize', 20)

% xlabel('Freq. (Hz)', 'FontSize', 16)

ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

%% Plotting zoomed spectra.

subplot(3, 1, 2)

boundedline((1:50)', WT_mean(1:50, :), prep_for_boundedline(WT_ci(1:50, :)))

hold on

plot([8 30; 8 30], repmat([all_dimensions(@min, WT_mean(1:50, :) - WT_ci(1:50, :)); all_dimensions(@max, WT_mean(1:50, :) + WT_ci(1:50, :))], 1, 2), '--r')

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

%% Plotting boxplots.

stats_name = dir([group_prefix, '_2.5_mins_ranksum_pct_power_*_stats.txt']);

broadband_stats = text_read(stats_name.name, '%s');

stats = nan(8, 12);

for b = 1:5
   
    band_stats_start = 24 + (b - 1)*13 + 1;
    
    band_stats_end = 23 + b*13;
    
    stats(b, :) = cellfun(@str2num, broadband_stats(band_stats_start:band_stats_end))';
    
end

stats_name = dir([group_prefix, '_8-30Hz_3-7cycles_6bands_2.5_mins_ranksum_pct_power_*_stats.txt']);

narrowband_stats = text_read(stats_name.name, '%s');

for b = 1:3
   
    band_stats_start = 24 + (b - 1)*13 + 1;
    
    band_stats_end = 23 + b*13;
    
    stats(5 + b, :) = cellfun(@str2num, narrowband_stats(band_stats_start:band_stats_end))';
    
end

median = stats(:, 5:6); Q1 = stats(:, 7:8); Q3 = stats(:, 9:10);

x = [(0:3:(3*7))' + 1 (0:3:(3*7))' + 2];

subplot(3,1,3)

errorbar(x, median, median - Q1, Q3 - median, 'd')

axis tight

y_limits = ylim;

y_range = diff(y_limits);

ylim([y_limits(1) - .1*y_range, y_limits(2) + .2*y_range])

xlim([0 25])

box off

ylabel({'Median & Quartiles'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)

set(gca, 'FontSize', 16, 'XTick', 1.5:3:23.5, 'XTickLabel', {'1-4Hz', '4-8Hz', '8-30Hz', '30-100Hz', '120-180Hz', '8-13Hz', '13-18Hz', '18-25Hz'})

save_as_eps(gcf, [group_prefix, '_spectra_boxplots'])

save_as_pdf(gcf, [group_prefix, '_spectra_boxplots'])

end