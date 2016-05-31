function make_all_band_figures

group_prefix = 'STR_M1';

short_chan_labels = {'Str.', 'M1'};

band_labels = {'1-4', '4-8', '8-12', '15-30', '40-100', '120-180'}; % , '0-200'};

no_bands = length(band_labels);

load([group_prefix, '_subjects.mat'])

figure

colorspec = [0 0 1; 0 .5 0];

%% Plotting period boxplots.
    
load([group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_pct_BP_high_2.5_min_secs_by_STR.mat'])

basetimes_mat = repmat(basetimes', [1, size(All_bp_max_start, 4), size(All_bp_max_start, 3)])/60;

bp_max_start = permute(All_bp_max_start(:, 1, :, :), [1 4 3 2])/(500*60) - basetimes_mat;

bp_max_end = permute(All_bp_max_end(:, 1, :, :), [1 4 3 2])/(500*60) - basetimes_mat;

for b = 1:no_bands
   
    subplot(no_bands, 6, (b - 1)*6 + 1)
    
    for fo = [1 2 4:length(folders)] % 1:length(folders)
        
        plot([-10 30], fo*[1 1], 'k')
        
        hold on
        
        for pd = 1:2
            
            plot([bp_max_start(fo, pd, b); bp_max_end(fo, pd, b)], fo*[1 1],... % '-d',...
                'LineWidth', 4, 'Color', colorspec(pd, :))
                % 'MarkerFaceColor', colorspec(pd, :), 'MarkerEdgeColor', colorspec(pd, :),...
                
            
        end
        
    end
    
    plot([0; 0], [1; length(folders)], ':k')
    
    xlim([-10 30]), ylim([.5 length(folders)+.5])
    
    box off
    
    if b == 1
        
        title({'Pd. Densest';'Str. BP'}, 'FontSize', 16)
        
    elseif b == no_bands
        
        xlabel({'Time'; '(m, Rel. Infusion)'}, 'FontSize', 16)
        
    end
    
    y = ylabel([band_labels{b}, ' Hz'], 'FontSize', 16, 'Rotation', 0);
    
    set(y, 'Units', 'Normalized', 'Position', [-0.7 0.4 0])
    
end

no_chans = length(chan_labels);

stat_suffixes = {'', 'pct_power_'}; no_stats = length(stat_suffixes);

stat_labels = {'Band Density', 'BP (% \Delta BL)'};

stats = nan(no_bands, 12, no_chans, no_stats);

%% Getting stats.

for st = 1:no_stats
    
    for ch = 1:no_chans
        
        % stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_t-test_', stat_suffixes{st}, chan_labels{ch}, '_stats.txt'];
        %
        % density_stats_cell = text_read(stats_name, '%s');
        
        for b = 1:no_bands
            
            stats_name = [group_prefix, '_1-200Hz_3-21cycles_7bands_kmeans_2.5_mins_', band_labels{b},...
                'Hz_individual_t-test_', stat_suffixes{st}, chan_labels{ch}, '_stats.txt'];
            
            density_stats_cell = text_read(stats_name, '%s');
            
            band_stats_start = 24 + length(folders)*13 + 1; % (b - 1)*13 + 1;
            
            band_stats_end = 23 + (length(folders) + 1)*13; % b*13;
            
            % band_stats_start = 24 + (b - 1)*13 + 1;
            % 
            % band_stats_end = 23 + b*13;

            stats(b, :, ch, st) = cellfun(@str2num, density_stats_cell(band_stats_start:band_stats_end))';
            
        end
        
    end
    
end

mean_stat = stats(:, 1:2, :, :); se_stat = stats(:, 3:4, :, :);

%% Plotting barplots.

for st = 1:no_stats
    
    for b = 1:no_bands
        
        subplot(no_bands, 6, (b - 1)*6 + st + 1)
        
        h = barwitherr(permute(se_stat(b, :, :, st), [3 2 1 4]), permute(mean_stat(b, :, :, st), [3 2 1 4]));
        
        set(h(1), 'FaceColor', [0 0 1]), set(h(2), 'FaceColor', [0 .5 0])
        
        bar_pos = get_bar_pos(h);
        
        p_vals = reshape(stats(b, 11, :, st), 2, 1);
        
        load('bonferroni_count')
        
        bar_pairs = {};
        
        for ch = 1:2
            
            if p_vals(ch) < .05 % /bonferroni_count
                
                bar_pairs = {bar_pairs{:}, [bar_pos(1, ch), bar_pos(2, ch)]};
                
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
        
        xlims = [bar_pos(1, 1) - 3*mean(diff(bar_pos)), bar_pos(2, 2) + 3*mean(diff(bar_pos))];
        
        xlim(xlims)
        
        set(gca, 'XTick', mean(bar_pos), 'XTickLabel', {'Str.', 'M1'})
        
        if b == 1
            
            title({stat_labels{st}; 'Mean \pm S.E.'}, 'FontSize', 16)
            
        elseif b == no_bands
            
            xlabel('Channel', 'FontSize', 16)
            
        end
        
    end
    
end

%% Plotting spectra.

group_prefixes = {'STR_w_M1', 'M1'};

for b = 1:no_bands
    
    for ch = 1:no_chans
        
        load([group_prefixes{ch}, '_1-200Hz_3-21cycles_7bands_kmeans_pct_', band_labels{b}, 'Hz_high_2.5_min_secs_pct_spectrum_no_130716_ch1_data_for_plot.mat'])
        
        All_mean_ci = norminv(1 - .05, 0, 1)*All_mean_se;
        
        [sig_lower, sig_higher] = find_sig(All_mean_mean, All_mean_ci);
        
        subplot(no_bands, 4, (b - 1)*4 + 2 + ch)
        
        boundedline((1:200)', All_mean_mean, prep_for_boundedline(All_mean_ci))
        
        axis tight
        
        add_stars(gca, (1:200)', logical(sig_lower), 0, [1 .5 0])
        
        add_stars(gca, (1:200)', logical(sig_higher), 1, [1 0 0])
        
        if b == 1
            
            title({[short_chan_labels{ch}, ' Pow. (% \Delta BL)']; 'Mean \pm 95% CI'}, 'FontSize', 16)
            
        elseif b == no_bands
        
            xlabel('Freq. (Hz)', 'FontSize', 16)
        
        end
        
    end
    
end

save_as_eps(gcf, [group_prefix, '_all_bands'])

save_as_pdf(gcf, [group_prefix, '_all_bands'])

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end