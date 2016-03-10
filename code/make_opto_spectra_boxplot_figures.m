function make_opto_spectra_boxplot_figures(group_prefix)

load([group_prefix, '_subjects.mat'], 'chan_labels', 'pd_labels')

no_chans = length(chan_labels);

no_pds_plotted = 2;

figure

%% Plotting broadband spectrum.

for ch = 1:no_chans
    
    load([group_prefix, '_pct_0-200Hz_10trials_pct_spectrum_ch', num2str(ch), '_data_for_plot.mat'])
    
    subplot(3, no_chans, ch)
    
    boundedline((1:200)', WT_mean(:, 1:no_pds_plotted), prep_for_boundedline(WT_ci(:, 1:no_pds_plotted)))
    
    axis tight
    
    set(gca, 'FontSize', 16)
    
    title(chan_labels{ch}, 'FontSize', 20)
    
    if ch == 1
        
        ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)
        
    else
        
        legend(pd_labels(1:no_pds_plotted))
        
    end
    
end

%% Plotting narrowband spectrum.
    
for ch = 1:no_chans

    load([group_prefix, '_8-30Hz_3-7cycles_6bands_8-30Hz_10trials_pct_spectrum_ch', num2str(ch), '_data_for_plot.mat'])
    
    subplot(3, no_chans, no_chans + ch)
    
    boundedline((8:.5:30)', WT_mean(:, 1:no_pds_plotted), prep_for_boundedline(WT_ci(:, 1:no_pds_plotted)))
    
    axis tight
    
    set(gca, 'FontSize', 16)
    
    % xlabel('Freq. (Hz)', 'FontSize', 16)
    
    if ch == 1
        
        ylabel({'Mean \pm 95% CI'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)
        
    end
    
end

%% Plotting boxplots.

for ch = 1:no_chans
    
    stats_name = dir([group_prefix, '_mins_ranksum_10trials_pct_power_', chan_labels{ch}, '_stats.txt']);
    
    if isempty(stats_name)
    
        stats_name = dir([group_prefix, '_mins_ranksum_10trials_power_pct_', chan_labels{ch}, '_stats.txt']);
        
    end
    
    broadband_stats = text_read(stats_name.name, '%s');
    
    stats = nan(8, 32);
    
    for b = 1:5
        
        band_stats_start = 54 + (b - 1)*33 + 1;
        
        band_stats_end = 53 + b*33;
        
        stats(b, :) = cellfun(@str2num, broadband_stats(band_stats_start:band_stats_end))';
        
    end
    
    stats_name = dir([group_prefix, '_8-30Hz_3-7cycles_6bands_mins_ranksum_10trials_pct_power_', chan_labels{ch}, '_stats.txt']);
    
    if isempty(stats_name)
    
        stats_name = dir([group_prefix, '_8-30Hz_3-7cycles_6bands_mins_ranksum_10trials_power_pct_', chan_labels{ch}, '_stats.txt']);
        
    end
    
    narrowband_stats = text_read(stats_name.name, '%s');
    
    for b = 1:3
        
        band_stats_start = 54 + (b - 1)*33 + 1;
        
        band_stats_end = 53 + b*33;
        
        stats(5 + b, :) = cellfun(@str2num, narrowband_stats(band_stats_start:band_stats_end))';
        
    end
    
    median = stats(:, 9:10); Q1 = stats(:, 13:14); Q3 = stats(:, 17:18);
    
    x = [(0:3:(3*7))' + 1 (0:3:(3*7))' + 2];
    
    subplot(3, no_chans, 2*no_chans + ch)
    
    errorbar(x, median, median - Q1, Q3 - median, 'd')
    
    axis tight
    
    y_limits = ylim;
    
    y_range = diff(y_limits);
    
    ylim([y_limits(1) - .1*y_range, y_limits(2) + .2*y_range])
    
    xlim([0 25])
    
    box off
    
    ylabel({'Median & Quartiles'; 'Power (% \Delta Baseline)'}, 'FontSize', 16)
    
    set(gca, 'FontSize', 16, 'XTick', 1.5:3:23.5, 'XTickLabel', {'1-4', '4-8', '8-30', '30-100', '120-180', '8-13', '13-18', '18-25'})
    
    xlabel('Freq. (Hz)', 'FontSize', 16)
    
end

save_as_eps(gcf, [group_prefix, '_spectra_boxplots'])

save_as_pdf(gcf, [group_prefix, '_spectra_boxplots'])

end