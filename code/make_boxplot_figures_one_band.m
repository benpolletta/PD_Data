function make_boxplot_figures_one_band(group_prefix, chan_label, band_index)

band_labels = {'1-4', '4-8', '8-15', '15-30', '40-100', '120-180'};

load([group_prefix, '_subjects.mat'], 'chan_labels')

figure

%% Plotting density boxplots.

stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_t-test_', chan_label, '_stats.txt'];

density_stats_cell = text_read(stats_name, '%s');

no_bands = 6;

density_stats = nan(no_bands, 12);

for b = 1:no_bands
   
    band_stats_start = 24 + (b - 1)*13 + 1;
    
    band_stats_end = 23 + b*13;
    
    density_stats(b, :) = cellfun(@str2num, density_stats_cell(band_stats_start:band_stats_end))';
    
end

mean_density = [density_stats(band_index, 1:2); nan nan];

se = [density_stats(band_index, 3:4); nan nan];

subplot(1, 2, 1)

h = barwitherr(se, mean_density);

set(h(1), 'FaceColor', [0 0 1]), set(h(2), 'FaceColor', [0 .5 0])

bar_pos = get_bar_pos(h);

p_vals = density_stats(:, 11);

load('bonferroni_count')
    
if p_vals(band_index) < .05/bonferroni_count
    
    bar_pairs = {[bar_pos(1, 1), bar_pos(2, 1)]};

    sigstar(bar_pairs, min(1, p_vals(band_index)*bonferroni_count))
    
end

axis tight

box off

y_lims = ylim; y_range = diff(y_lims);

if y_lims(1) ~= 0

    ylim([y_lims(1) - .1*y_range, y_lims(2) + .1*y_range])
    
else
    
    ylim([y_lims(1), y_lims(2) + .1*y_range])
    
end

xlims = [bar_pos(1, 1) - 1.5*diff(bar_pos(:, 1)), bar_pos(2, 1) + 1.5*diff(bar_pos(:, 1))];

xlim(xlims)

% ylabel({'Density (Mean \pm S.E.)'}, 'FontSize', 16)

set(gca, 'XTick', bar_pos(:, 1)', 'XTickLabel', {'Pre', 'Post'}, 'FontSize', 16)

xticklabel_rotate([], 45, [], 'FontSize', 16)

%% Plotting power barplots.

stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_t-test_pct_power_', chan_label, '_stats.txt'];

power_stats_cell = text_read(stats_name, '%s');

for b = 1:no_bands
   
    band_stats_start = 24 + (b - 1)*13 + 1;
    
    band_stats_end = 23 + b*13;
    
    power_stats(b, :) = cellfun(@str2num, power_stats_cell(band_stats_start:band_stats_end))';
    
end

mean_power = [power_stats(band_index, 1:2); nan nan];

se = [power_stats(band_index, 3:4); nan nan];

subplot(1, 2, 2)

h = barwitherr(se, mean_power);

set(h(1), 'FaceColor', [0 0 1]), set(h(2), 'FaceColor', [0 .5 0])

bar_pos = get_bar_pos(h);

p_vals = power_stats(:, 11);

load('bonferroni_count')
    
if p_vals(band_index) < .05/bonferroni_count
    
    bar_pairs = {[bar_pos(1, 1), bar_pos(2, 1)]};
    
    sigstar(bar_pairs, min(1, p_vals(band_index)*bonferroni_count))
    
end
    
axis tight

box off

y_lims = ylim; y_range = diff(y_lims);

if y_lims(1) ~= 0

    ylim([y_lims(1) - .1*y_range, y_lims(2) + .1*y_range])
    
else
    
    ylim([y_lims(1), y_lims(2) + .1*y_range])
    
end

xlims = [bar_pos(1, 1) - 1.5*diff(bar_pos(:, 1)), bar_pos(2, 1) + 1.5*diff(bar_pos(:, 1))];

xlim(xlims)

% ylabel({'Power (%\Delta BL, Mean \pm S.E.)'}, 'FontSize', 16)

set(gca, 'XTick', bar_pos(:, 1)', 'XTickLabel', {'Pre', 'Post'}, 'FontSize', 16)

xticklabel_rotate([], 45, [], 'FontSize', 16)

% xlabel('Freq. Range (Hz)', 'FontSize', 16)

save_as_eps(gcf, [group_prefix, '_', chan_label, '_', band_labels{band_index}, '_boxplots'])

save_as_pdf(gcf, [group_prefix, '_', chan_label, '_', band_labels{band_index}, '_boxplots'])

end