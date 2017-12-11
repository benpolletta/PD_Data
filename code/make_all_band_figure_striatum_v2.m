function make_all_band_figure_striatum_v2

peak_suffix = '_kmeans_win_420_1020';

channel_prefix = 'STR_w_M1';

channel_title = 'Striatal';

freq_limit = 100;

p_val = .05;

test_flag = 'ttest';

short_group_labels = {'All', 'M1+', 'M1-'};

band_number_flag =  {'7', '8', '8', '7', '7', '7'}; % {'7', '7', '7', '7', '7', '7'}; %

band_labels = {'1-4', '5-8', '9-14', '15-30', '40-100', '120-180'}; % , '0-200'}; % '4-8', '8-12', '15-30', '40-100', '120-180'}; % 

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

colorspec = [0 0 1; 0 .5 0];

groups = {'missing_2'}; %, 'M1_not_increased', 'M1_increased'}; % , 'M1_increased', 'M1_not_increased'};

no_groups = length(groups);

folder_chi = ones(length(folders), length(groups)); % Indicator function for which folders are in each group.

for g = 1:length(groups)
    
    folder_cell = load(groups{g}); % Getting cell of excluded filenames for each group.
    
    folder_cell = getfield(folder_cell, groups{g});
    
    folder_cell = folder_cell{2};
    
    for fo = 1:length(folder_cell)
        
        folder_chi(strcmp(folders, folder_cell{fo}), g) = 0; % Setting folder_chi to zero for each excluded folder.
        
    end
    
end

%% Plotting spectra.

figure

freq_lims = [1 50; 50 200];
frequencies = 1:200;
no_fls = size(freq_lims, 1);

tight_subplot_handles = tight_subplot(no_bands, no_fls);

for b = 1:no_bands
    
    % Plotting without stars.
    
    for g = 1:no_groups  % Loading & plotting mean & SE data for each group.
        
        for fl = 1:no_fls
            
            freq_index = frequencies >= freq_lims(fl, 1) & frequencies <= freq_lims(fl, 2);
            
            load([channel_prefix, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_', groups{g}, '_ch1_data_for_plot.mat'])
            
            All_mean_ci = norminv(1 - max(p_val), 0, 1)*All_mean_se; % Constructing CI from p_val & SE.
            
            axes(tight_subplot_handles((b - 1)*no_groups*no_fls + (g - 1)*no_fls + fl))
            
            boundedline(frequencies(freq_index)', All_mean_mean(freq_index, :), prep_for_boundedline(All_mean_ci(freq_index, :)))
            
            axis tight
            
            y_lims((g - 1)*no_fls + fl, :) = ylim; % Collecting y limits to make them uniform across plots.
            
        end
        
    end
    
    y_extremes = [min(y_lims(:, 1)) max(y_lims(:, 2))];
    
    load([channel_prefix, '_1-200Hz_3-21cycles_', band_number_flag{b}, 'bands', peak_suffix, '_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat']) % Loading mean data for all individuals.
    
    for g = 1:no_groups
        
        for fl = 1:no_fls
            
            freq_index = frequencies >= freq_lims(fl, 1) & frequencies <= freq_lims(fl, 2);
            
            p_vals = nan(sum(freq_index), 2);
            
            for f = frequencies(freq_index)
                
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
            
            axes(tight_subplot_handles((b - 1)*no_groups*no_fls + (g - 1)*no_fls + fl))
            
            set(gca, 'ylim', y_extremes)
            
            add_stars_one_line(gca, frequencies(freq_index)', logical(test(freq_index, :, 1)), 0, 'c_order', colors(:, :, 1)) % logical(test(:, 1)), 0, [1 .5 0])
            
            add_stars_one_line(gca, frequencies(freq_index)', logical(test(freq_index, :, 2)), 1, 'c_order', colors(:, :, 2)) % logical(test(:, 2)), 1, [1 0 0])
            
            if b == 1 && g == 1 && fl == 1
                
                legend({'Pre-Infusion', 'Post-Infusion'})
                
                if length(groups) > 1
                    
                    title({[short_group_labels{g}, ' ', channel_title, ' Pow. (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10)
                    
                else
                    
                    title({[channel_title, ' Pow. (% \Delta BL)']; ['Mean \pm ', num2str(100*(1 - 2*max(p_val)), '%g'), '% CI']}, 'FontSize', 10)
                    
                end
                
            end
            
            set(gca, 'FontSize', 10, 'XTickLabelMode', 'auto', 'YTickLabelMode', 'auto')
                
            if b == no_bands
                
                xlabel('Freq. (Hz)', 'FontSize', 10)
                
            end
            
            if g == 1 && fl == 1
            
                y = ylabel({band_labels{b}; 'Hz BPD'}, 'FontSize', 10, 'Rotation', 45);
                
            end
            
            if fl == 1
                
                set(gca, 'XTick', [15 30 45])
                
            elseif fl > 1
                
                set(gca, 'YTick', [], 'YColor', 'w')
                
            end
            
        end
    
        sync_axes(tight_subplot_handles((b - 1)*no_groups*no_fls + (g - 1)*no_fls + (1:2)), 'y')
        
    end
    
end

name = ['striatum', peak_suffix, '_all_bands_spec_allfreqs', p_tag, '_', test_flag];

set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', .9*[0 0 9.1 18.2], 'PaperUnits', 'centimeters', 'PaperPosition', .9*[0 0 9.1 18.2])

print(gcf, '-painters', '-dpdf', '-r600', [name, '.pdf'])

print(gcf, '-painters', '-depsc', '-r600', [name, '.eps'])

saveas(gcf, [name, '.fig'])

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