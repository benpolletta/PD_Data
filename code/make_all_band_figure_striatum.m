function make_all_band_figure_striatum

channel_prefix = 'STR_w_M1';

channel_title = 'Striatal';

freq_limit = 100;

p_val = .02/2;

test_flag = 'ttest';

short_group_labels = {'All', 'M1+', 'M1-'};

band_number_flag = {'7', '8', '8', '7', '7', '7'};

band_labels = {'1-4', '5-8', '9-14', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

colorspec = [0 0 1; 0 .5 0];

group = {'missing_2'}; %, 'M1_not_increased', 'M1_increased'}; % , 'M1_increased', 'M1_not_increased'};

no_groups = length(group);

folder_chi = ones(length(folders), length(group)); % Indicator function for which folders are in each group.

for g = 1:length(group)
    
    folder_cell = load(group{g}); % Getting cell of excluded filenames for each group.
    
    folder_cell = getfield(folder_cell, group{g});
    
    folder_cell = folder_cell{2};
    
    for fo = 1:length(folder_cell)
        
        folder_chi(strcmp(folders, folder_cell{fo}), g) = 0; % Setting folder_chi to zero for each excluded folder.
        
    end
    
end

%% Plotting spectra.

for b = 1:no_bands
    
    % Computing stats between post-infusion increases.

    % load('missing_2.mat')
    % 
    % missing_2 = missing_2{2};
    % 
    % no_excluded = length(missing_2);
    % 
    % folder_index = ones(1, length(folders));
    % 
    % for e = 1:no_excluded
    % 
    %     folder_index = folder_index - strcmp(folders, missing_2{e});
    % 
    % end
    
    % Plotting without stars.
    
    for g = 1:no_groups  % Loading & plotting mean & SE data for each group.
        
        load([channel_prefix, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_', group{g}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se; % Constructing CI from p_val & SE.
        
        subplot(no_bands, 1, b)
        
        boundedline((1:freq_limit)', All_mean_mean(1:freq_limit, :), prep_for_boundedline(All_mean_ci(1:freq_limit, :)))
        
        axis tight
        
        y_lims(g, :) = ylim; % Collecting y limits to make them uniform across plots.
        
    end
    
    y_extremes = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    load([channel_prefix, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat']) % Loading mean data for all individuals.
    
    for g = 1:no_groups
        
        p_vals = nan(freq_limit, 2);
        
        for f = 1:freq_limit
            
            if strcmp(test_flag, 'ranksum')
                
                p_vals(f, 1) = ranksum(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'right');
                
                p_vals(f, 2) = ranksum(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'left');
                
            elseif strcmp(test_flag, 'ttest')
                
                [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'left');
                
            elseif strcmp(test_flag, 'ttest2')
                
                [~, p_vals(f, 1)] = ttest2(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest2(All_mean(f, logical(folder_chi(:, g)), 1)', All_mean(f, logical(folder_chi(:, g)), 2)', 'tail', 'left');
                
            end
            
        end
        
        [test, colors, p_tag] = test_p_vals(p_vals, p_val, [1 .5 0; 1 0 0]); % Getting test matrix, gradient of colors if p_val is a vector.
        
        % test = p_vals < p_val;
        
        h = subplot(no_bands, 1, b);
        
        set(h, 'ylim', y_extremes)
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 1)), 0, colors(:, :, 1)) % logical(test(:, 1)), 0, [1 .5 0])
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, :, 2)), 1, colors(:, :, 2)) % logical(test(:, 2)), 1, [1 0 0])
        
        % % Finding nonoverlap of confidence intervals.
        % load([channel_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_', group{g}, '_ch1_data_for_plot.mat'])
        % 
        % All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se;
        % 
        % [sig_lower, sig_higher] = find_sig(All_mean_mean(1:freq_limit, :), All_mean_ci(1:freq_limit, :));
        % 
        % add_stars(gca, (1:freq_limit)', logical(sig_lower), 0, [1 .5 0])
        % 
        % add_stars(gca, (1:freq_limit)', logical(sig_higher), 1, [1 0 0])
        
        if b == 1
            
            legend({'Pre-Infusion', 'Post-Infusion'})
            
            if length(group) > 1
            
                title({[short_group_labels{g}, ' ', channel_title, ' Pow. (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10)
                
            else
            
                title({[channel_title, ' Pow. (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10)
            
            end
                
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 10)
        
        end
        
        y = ylabel({band_labels{b}; 'Hz BPD'}, 'FontSize', 10, 'Rotation', 45); % , 'Rotation', 0);
        
        % set(y, 'Units', 'Normalized', 'Position', [-0.1 0.5 0])
        
    end
    
end

set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', [0 0 9.1 1.5*9.1], 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 9.1 1.5*9.1])

print(gcf, '-painters', '-dpdf', '-r600', ['striatum_all_bands_spec_f', num2str(freq_limit), p_tag, '_', test_flag, '.pdf'])

print(gcf, '-painters', '-depsc', '-r600', ['striatum_all_bands_spec_f', num2str(freq_limit), p_tag, '_', test_flag, '.eps'])

saveas(gcf, ['striatum_all_bands_spec_f', num2str(freq_limit), p_tag, '_', test_flag, '.fig'])

end

function [test, colors, p_tag] = test_p_vals(p_vals, p_val, colors_in)

no_ps = length(p_val);

no_tests = size(p_vals, 2);

test = nan([size(p_vals, 1), no_ps, no_tests]);

colors = nan(no_ps, 3, no_tests);
    
for t = 1:no_tests
    
    test(:, 1, t) = p_vals(:, t) < p_val(1);
    
    for p = 2:no_ps
    
        test(:, p, t) = p_vals(:, t) >= p_val(p - 1) & p_vals(:, t) < p_val(p);
    
    end

    colors(:, :, t) = flipud(color_gradient(no_ps, .5*colors_in(t, :), colors_in(t, :)));
    
end

if isscalar(p_val)
    
    p_tag = sprintf('_p%g', p_val);
    
else
    
    p_tag = sprintf('_p%gto%g', p_val(1), p_val(end));

end

end