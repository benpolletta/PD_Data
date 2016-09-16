function make_supplemental_figure

p_val = .01; freq_limit = 100;

short_chan_labels = {'M1'};

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

colorspec = [0 0 1; 0 .5 0];

load('missing_2')

%% Plotting period boxplots.
    
load('STR_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_BP_high_2.5_min_secs.mat') % _by_STR.mat'])

basetimes_mat = repmat(basetimes', [1, size(All_bp_max_start, 4), size(All_bp_max_start, 3)])/60;

[bp_max_start, bp_max_end] = deal(nan(length(folders), 2, no_bands + 1));

for s_id = 1:2
    
    s_ind = striatal_id == s_id;
    
    bp_max_start(s_ind, :, :) = permute(All_bp_max_start(s_ind, s_id, :, :), [1 4 3 2])/(500*60) - basetimes_mat(s_ind, :, :);

    bp_max_end(s_ind, :, :) = permute(All_bp_max_end(s_ind, s_id, :, :), [1 4 3 2])/(500*60) - basetimes_mat(s_ind, :, :);
    
end

folder_chi = ones(size(folders));

folder_cell = missing_2{2};

for fo = 1:length(folder_cell)
    
    folder_chi(strcmp(folders, folder_cell{fo})) = 0;
    
end

M1_beta = ones(size(folders));

load('M1_groups')

for fo = 1:length(M1_groups{2})
   
    M1_beta(strcmp(folders, M1_groups{2}{fo})) = 0;
    
end

subj_index = find(folder_chi);

mean_bp_max_start = mean(bp_max_start(subj_index, :, :));

std_bp_max_start = std(bp_max_start(subj_index, :, :)); 

for b = 1:no_bands

    fprintf('Mean Start of Max. High %s Hz Density = %f\n', band_labels{b}, mean_bp_max_start(:, 2, b))

    fprintf('St. Dev. Start of Max. High %s Hz Density = %f\n', band_labels{b}, std_bp_max_start(:, 2, b))
   
    subplot(no_bands, 3, (b - 1)*3 + 1)
    
    folder_index = 1;
    
    for fo = subj_index % 1:length(folders)
        
        plot([-10 30], folder_index*[1 1], 'k')
        
        hold on
        
        for pd = 1:2
            
            if M1_beta(fo)
                
                plot([bp_max_start(fo, pd, b); bp_max_end(fo, pd, b)], folder_index*[1 1],... % '-d',...
                    'LineWidth', 4.5, 'Color', colorspec(pd, :))
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            else
                
                plot([bp_max_start(fo, pd, b); bp_max_end(fo, pd, b)], folder_index*[1 1],... % '-d',...
                    'LineWidth', 2.5, 'Color', colorspec(pd, :)) % , 'LineStyle', '--')
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            end
            
        end
        
        folder_index = folder_index + 1;
        
    end
    
    plot([0; 0], [1; length(folders)], 'k')
    
    xlim([-10 30]), ylim([.5 (length(subj_index) +.5)])
    
    box off
    
    if b == 1
        
        title({'Pd. of Highest';'Striatal BPD'}, 'FontSize', 10)
        
    elseif b == no_bands
        
        xlabel({'Time'; '(m, Rel. Infusion)'}, 'FontSize', 10)
        
    end
    
    y = ylabel({band_labels{b}; 'Hz'}, 'FontSize', 10, 'Rotation', 45);
    
    % set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
    
end

no_chans = length(chan_labels);

%% Plotting spectra.

chan_prefixes = {'M1'}; no_chans = length(chan_prefixes);

for b = 1:no_bands
    
    for ch = 1:no_chans
        
        load([chan_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_missing_2_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
        
        subplot(no_bands, 2 + no_chans, (b - 1)*(2 + no_chans) + 1 + ch)
        
        boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
        
        axis tight
        
        load([chan_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_data_for_plot.mat'])
        
        p_vals = nan(freq_limit, 2);
        
        for f = 1:freq_limit
                
                [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'right');
                
                [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'left');
            
        end
        
        test = p_vals < p_val;
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 1)), 0, [1 .5 0])
        
        add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 2)), 1, [1 0 0])
        
        if b == 1
            
            title(short_chan_labels{ch}, 'FontSize', 10)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 10)
        
        end
        
    end
    
end

%% Plotting PLV.

for b = 1:no_bands
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_missing_2_Coh_sec_pct_data_for_plot.mat'])
    
    All_mean_ci = norminv(1 - p_val, 0, 1)*All_mean_se;
    
    subplot(no_bands, 2 + no_chans, (b - 1)*(2 + no_chans) + 2 + no_chans)
    
    boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
    
    axis tight
    
    load(['STR_M1_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV_data_for_plot.mat'])
    
    p_vals = nan(freq_limit, 2);
    
    for f = 1:freq_limit
        
        [~, p_vals(f, 1)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'right');
        
        [~, p_vals(f, 2)] = ttest(All_mean(f, logical(folder_chi), 1)', All_mean(f, logical(folder_chi), 2)', 'tail', 'left');
        
    end
    
    test = p_vals < p_val;
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 1)), 0, [1 .5 0])
    
    add_stars_one_line(gca, (1:freq_limit)', logical(test(:, 2)), 1, [1 0 0])
    
    if b == 1
        
        title('PLV', 'FontSize', 10)
        
    elseif b == no_bands
        
        xlabel('Freq. (Hz)', 'FontSize', 10)
        
    end
    
end
    
set(gcf, 'PaperOrientation', 'landscape', 'Units', 'centimeters', 'Position', [0 0 18.2 9.1*1.5], 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 18.2 9.1*1.5])

print(gcf, '-painters', '-dpdf', '-r600', 'supplementary_figure.pdf')

print(gcf, '-painters', '-depsc', '-r600', 'supplementary_figure.eps')
    
end