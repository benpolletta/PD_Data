function make_boxplot_figures(group_prefix, chan_label)

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'};

load([group_prefix, '_subjects.mat'], 'chan_labels', 'folders')

figure

%% Plotting density boxplots.

% stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_t-test_', chan_label, '_stats.txt'];
%
% density_stats_cell = text_read(stats_name, '%s');

no_bands = 6;

[density_stats, power_stats] = deal(nan(no_bands, 12));

for b = 1:no_bands

    stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_', band_labels{b}, 'Hz_individual_t-test_', chan_label, '_stats.txt'];
   
    density_stats_cell = text_read(stats_name, '%s');
    
    band_stats_start = 24 + length(folders)*13 + 1; % (b - 1)*13 + 1;
    
    band_stats_end = 23 + (length(folders) + 1)*13; % b*13;
    
    density_stats(b, :) = cellfun(@str2num, density_stats_cell(band_stats_start:band_stats_end))';
    
end

mean_density = density_stats(:, 1:2); se = density_stats(:, 3:4);

x = [(0:3:(3*5))' + 1 (0:3:(3*5))' + 2];

subplot(2,1,1)

h = barwitherr(se, mean_density);

set(h(1), 'FaceColor', [0 0 1]), set(h(2), 'FaceColor', [0 .5 0])

bar_pos = get_bar_pos(h);

p_vals = density_stats(:, 11);

load('bonferroni_count')

bar_pairs = {};

for b = 1:no_bands
    
    if p_vals(b) < .05 % /bonferroni_count
        
        bar_pairs = {bar_pairs{:}, [bar_pos(1, b), bar_pos(2, b)]};
        
    end
    
end

sigstar(bar_pairs, p_vals(p_vals < .05)) % /bonferroni_count))

axis tight

box off

y_lims = ylim; y_range = diff(y_lims);

if y_lims(1) ~= 0

    ylim([y_lims(1) - .1*y_range, y_lims(2) + .1*y_range])
    
else
    
    ylim([y_lims(1), y_lims(2) + .1*y_range])
    
end

xlims = [bar_pos(1, 1) - 3*mean(diff(bar_pos)), bar_pos(2, no_bands) + 3*mean(diff(bar_pos))];

xlim(xlims)

ylabel({'Mean \pm S.E.'; 'Density'}, 'FontSize', 16)

set(gca, 'XTick', mean(bar_pos), 'XTickLabel', band_labels, 'FontSize', 16)

xticklabel_rotate([], 45, [], 'FontSize', 16)

%% Plotting power barplots.

% stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_t-test_pct_power_', chan_label, '_stats.txt'];

% power_stats_cell = text_read(stats_name, '%s');

for b = 1:no_bands

    stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_', band_labels{b}, 'Hz_individual_t-test_pct_power_', chan_label, '_stats.txt'];
   
    power_stats_cell = text_read(stats_name, '%s');
    
    band_stats_start = 24 + length(folders)*13 + 1; % (b - 1)*13 + 1;
    
    band_stats_end = 23 + (length(folders) + 1)*13; % b*13;
    
    power_stats(b, :) = cellfun(@str2num, power_stats_cell(band_stats_start:band_stats_end))';
    
end

mean_power = power_stats(:, 1:2); se = power_stats(:, 3:4);

x = [(0:3:(3*5))' + 1 (0:3:(3*5))' + 2];

subplot(2,1,2)

h = barwitherr(se, mean_power);

set(h(1), 'FaceColor', [0 0 1]), set(h(2), 'FaceColor', [0 .5 0])

bar_pos = get_bar_pos(h);

p_vals = power_stats(:, 11);

load('bonferroni_count')

bar_pairs = {};

for b = 1:no_bands
    
    if p_vals(b) < .05 % /bonferroni_count
        
        bar_pairs = {bar_pairs{:}, [bar_pos(1, b), bar_pos(2, b)]};
        
    end
    
end

sigstar(bar_pairs, p_vals(p_vals < .05)) % /bonferroni_count))

axis tight

box off

y_lims = ylim; y_range = diff(y_lims);

if y_lims(1) ~= 0

    ylim([y_lims(1) - .1*y_range, y_lims(2) + .1*y_range])
    
else
    
    ylim([y_lims(1), y_lims(2) + .1*y_range])
    
end

xlims = [bar_pos(1, 1) - 3*mean(diff(bar_pos)), bar_pos(2, no_bands) + 3*mean(diff(bar_pos))];

xlim(xlims)

ylabel({'Mean \pm S.E.'; 'Power (%\Delta BL)'}, 'FontSize', 16)

set(gca, 'XTick', mean(bar_pos), 'XTickLabel', band_labels, 'FontSize', 16)

xticklabel_rotate([], 45, [], 'FontSize', 16)

xlabel('Freq. Range (Hz)', 'FontSize', 16)

save_as_eps(gcf, [group_prefix, '_', chan_label, '_boxplots'])

save_as_pdf(gcf, [group_prefix, '_', chan_label, '_boxplots'])

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end