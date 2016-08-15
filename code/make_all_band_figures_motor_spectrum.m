function make_all_band_figures_motor_spectrum

channel_prefix = 'M1';

short_group_labels = {'All M1,', 'M1 \beta\leq,', 'M1 \beta>,'};

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load('STR_M1_subjects.mat')

figure

colorspec = [0 0 1; 0 .5 0];

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

folder_cell = {'130716', '130830'};

for fo = 1:length(folder_cell)
    
    folder_chi(strcmp(folders, folder_cell{fo})) = 0;
    
end

subj_index = find(folder_chi);

mean_bp_max_start = mean(bp_max_start(subj_index, :, :));

std_bp_max_start = std(bp_max_start(subj_index, :, :)); 

for b = 1:no_bands

    fprintf('Mean Start of Max. High %s Hz Density = %f\n', band_labels{b}, mean_bp_max_start(:, 2, b))

    fprintf('St. Dev. Start of Max. High %s Hz Density = %f\n', band_labels{b}, std_bp_max_start(:, 2, b))
   
    subplot(no_bands, 4, (b - 1)*4 + 1)
    
    folder_index = 1;
    
    for fo = subj_index % 1:length(folders)
        
        plot([-10 30], folder_index*[1 1], 'k')
        
        hold on
        
        for pd = 1:2
            
            plot([bp_max_start(fo, pd, b); bp_max_end(fo, pd, b)], folder_index*[1 1],... % '-d',...
                'LineWidth', 4, 'Color', colorspec(pd, :))
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            
        end
        
        folder_index = folder_index + 1;
        
    end
    
    plot([0; 0], [1; length(folders)], ':k')
    
    xlim([-10 30]), ylim([.5 (length(subj_index) +.5)])
    
    box off
    
    if b == 1
        
        title({'Pd. Densest';'Str. BP'}, 'FontSize', 16)
        
    elseif b == no_bands
        
        xlabel({'Time'; '(m, Rel. Infusion)'}, 'FontSize', 16)
        
    end
    
    y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
    
    set(y, 'Units', 'Normalized', 'Position', [-0.35 0.4 0])
    
end

%% Plotting spectra.

group = {'missing_2', 'M1_increased', 'M1_not_increased'};

no_groups = length(group);

for b = 1:no_bands
    
    for g = 1:no_groups
        
        load(['M1_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_', group{g}, '_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;
        
        [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        
        subplot(no_bands, 4, (b - 1)*4 + 1 + g)
        
        boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
        
        axis tight
        
        add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
        
        add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
        
        if b == 1
            
            title({[short_group_labels{g}, ' Pow. (% \Delta BL)']; 'Mean \pm 95% CI'}, 'FontSize', 16)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 16)
        
        end
        
    end
    
end

% %% Plotting PLV.
% 
% for b = 1:no_bands
%     
%     load([group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_PLV', missing{1}, '_Coh_sec_pct_data_for_plot.mat'])
%     
%     % All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;
%     
%     [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
%     
%     subplot(no_bands, 4, (b - 1)*4 + 4)
%     
%     boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
%     
%     axis tight
%     
%     add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
%     
%     add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
%     
%     if b == 1
%         
%         title({'Phase-Locking Mag. (% \Delta BL)'; 'Mean \pm 95% CI'}, 'FontSize', 16)
%         
%     elseif b == no_bands
%         
%         xlabel('Freq. (Hz)', 'FontSize', 16)
%         
%     end
%     
% end

save_as_eps(gcf, [channel_prefix, '_all_bands_motor_spec'])

save_as_pdf(gcf, [channel_prefix, '_all_bands_motor_spec'])

end