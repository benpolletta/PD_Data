function make_all_band_figures_transmission_by_groups(p_val)

channel_prefixes = {'STR_w_M1', 'M1'};

channel_labels = {'Striatum', 'M1'};

short_group_labels = {'M1 \beta\leq', 'M1 \beta>'};

group_flags = {'M1_increased', 'M1_not_increased'}; no_groups = length(group_flags);

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

%% Computing stats between channels (post-infusion increases).

load('M1_groups.mat')

folder_index = cell(1, 2);

for g = 1:2
    
    no_excluded = length(M1_groups{g});
    
    folder_index{g} = ones(1, length(folders));
    
    for e = 1:no_excluded
    
        folder_index{g} = folder_index{g} - strcmp(folders, M1_groups{g}{e});
        
    end

end

freq_limit = 200;

for b = 1:no_bands
    
    All_mean_channels = nan(200, length(folder_index{1}), 2);
    
    for ch = 1:2
        
        load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
        
        All_mean_channels(:, :, ch) = All_mean(:, :, 2);
        
    end
    
    for g = 1:2
        
        mean_mat = nan(freq_limit, 2);
        
        for ch = 1:2
            
            load([channel_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_', group_flags{g}, '_ch1_data_for_plot.mat'])
            
            mean_mat(:, ch) = All_mean_mean(:, 2);
            
        end
        
        subplot(no_bands, no_groups, (b - 1)*no_groups + g)
        
        plot((1:freq_limit)', -diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)
        
        axis tight
        
        y_lims(g, :) = ylim;
        
    end
    
    y_extremes = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    for g = 1:2
        
        p_vals = nan(freq_limit, 2);
        
        for f = 1:freq_limit
            
            [~, p_vals(f, 1)] = ttest2(All_mean_channels(f, logical(folder_index{g}), 1)', All_mean(f, logical(folder_index{g}), 2)', 'tail', 'left');
            
            [~, p_vals(f, 2)] = ttest2(All_mean_channels(f, logical(folder_index{g}), 1)', All_mean(f, logical(folder_index{g}), 2)', 'tail', 'right');
            
        end
        
        test = p_vals < p_val;
        
        h = subplot(no_bands, no_groups, (b - 1)*no_groups + g);
        
        set(h, 'ylim', y_extremes)
        
        add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [1 1 0])
        
        add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [.5 0 .5])
        
        if ~any(test ~= 0), hold on, end
        
        plot((1:200)', zeros(200, 1), 'k:')
        
        if b == 1
            
            title({[short_group_labels{g}, ' Post-Inf. Pow. (%\Delta BL)']; [channel_labels{1}, ' - ', channel_labels{2}]}, 'FontSize', 16);
            
        end
        
        if g == 1
            
            y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
            
            set(y, 'Units', 'Normalized', 'Position', [-0.2 0.4 0])
            
        end
        
        box off
        
    end
    
end

save_as_eps(gcf, sprintf('STR_M1_all_bands_spec_channel_comp_p%g', p_val))

save_as_pdf(gcf, sprintf('STR_M1_all_bands_spec_channel_comp_p%g', p_val))

end