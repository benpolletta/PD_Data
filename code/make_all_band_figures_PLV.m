function make_all_band_figures_PLV(freq_limit, p_val)

channel_prefix = 'STR_M1';

short_group_labels = {'M1 \beta\leq', 'M1 \beta>'};

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

colorspec = [0 0 1; 0 .5 0];

%% Plotting spectra.

group = {'M1_increased', 'M1_not_increased'};

no_groups = length(group);

for b = 1:no_bands
    
    for g = 1:no_groups
        
        load([channel_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_', group{g}, '_Coh_sec_pct_data_for_plot.mat'])
        
        subplot(no_bands, no_groups + 1, (b - 1)*(no_groups + 1) + g)
        
        boundedline((1:freq_limit)', All_mean_mean(1:freq_limit, :), prep_for_boundedline(All_mean_ci(1:freq_limit, :)))
        
        axis tight
        
        y_lims(g, :) = ylim;
        
    end
    
    y_extremes = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    for g = 1:no_groups
        
        h = subplot(no_bands, no_groups + 1, (b - 1)*(no_groups + 1) + g);
        
        set(h, 'ylim', y_extremes)
        
        load([channel_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_', group{g}, '_Coh_sec_pct_data_for_plot.mat'])
        
        [sig_lower, sig_higher] = find_sig(All_mean_mean(1:freq_limit, :), All_mean_ci(1:freq_limit, :));
        
        add_stars(gca, (1:freq_limit)', logical(sig_lower), 0, [1 .5 0])
        
        add_stars(gca, (1:freq_limit)', logical(sig_higher), 1, [1 0 0])
        
        if b == 1
            
            title({['PLV (% \Delta BL), ', short_group_labels{g}]; 'Mean \pm 95% CI'}, 'FontSize', 16)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 16)
        
        end
        
        if g == 1
            
            y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
            
            set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
            
        end
        
        box off
        
    end
    
end

%% Computing stats between post-infusion increases.

load('M1_groups.mat')

folder_index = cell(1, 2);

for g = 1:2
    
    no_excluded = length(M1_groups{g});
    
    folder_index{g} = ones(1, length(folders));
    
    for e = 1:no_excluded
    
        folder_index{g} = folder_index{g} - strcmp(folders, M1_groups{g}{e});
        
    end

end

for b = 1:no_bands
    
    load([channel_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit
        
        [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_index{2}), 2)', All_mean(f, logical(folder_index{1}), 2)', 'tail', 'left');
        
        [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_index{2}), 2)', All_mean(f, logical(folder_index{1}), 2)', 'tail', 'right');
        
    end
    
    test = p_vals < p_val;
    
    mean_mat = nan(freq_limit, 2);
    
    for g = 1:2
        
        load([channel_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_', group{g}, '_Coh_sec_pct_data_for_plot.mat'])
        
        mean_mat(:, g) = All_mean_mean(1:freq_limit, 2);
        
    end
        
    subplot(no_bands, no_groups + 1, b*(no_groups + 1))
    
    plot((1:freq_limit)', diff(mean_mat(1:freq_limit, :), [], 2), 'k', 'LineWidth', 1.5)
    
    axis tight
    
    add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 1)), 0, [0 1 1])
    
    add_stars(gca, (1:freq_limit)', logical(test(1:freq_limit, 2)), 1, [1 0 1])
    
    if ~any(test ~= 0), hold on, end
    
    plot((1:freq_limit)', zeros(freq_limit, 1), 'k:')
    
    if b == 1
        
        title({['Post-Inf. PLV (%\Delta BL)']; [short_group_labels{2}, ' - ', short_group_labels{1}]}, 'FontSize', 16);
        
    end
    
    box off
    
end

save_as_eps(gcf, [channel_prefix, '_all_bands_PLV_f', num2str(freq_limit), '_p', num2str(p_val, '%g')])

save_as_pdf(gcf, [channel_prefix, '_all_bands_PLV_f', num2str(freq_limit), '_p', num2str(p_val, '%g')])

end