function PD_beta_blocks_rel_infusion_laser_plots(subject_mat, measure, no_trials)
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

pd_indices = cell(no_folders, 1);

if isempty(no_trials)
    
    trials = cell(no_folders, no_pds);
    
    no_trials = nan(no_folders, no_pds);
    
    for fo = 1:no_folders
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        load([subj_name, '_wav_laser_artifacts.mat'])
        
        pd_indices{fo} = logical(laser_periods);
        
        for pd = 1:no_pds
            
            trials{fo, pd} = index_to_blocks(laser_periods(:, pd));
            
            no_trials(fo, pd) = size(trials, 1);
            
        end
        
    end
    
    max_no_trials = all_dimensions(@max, no_trials);
    
else
    
    max_no_trials = no_trials;
    
end

no_secs = max_no_trials*5;
    
load([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_', num2str(no_trials), 'trials', measure, '.mat'])

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

p_vals = nan(no_comparisons, no_bands, no_chans, 2);

for b = 1:no_bands
    
    figure
    
    for ch = 1:no_chans
        
        pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
        
        for pd = 1:no_pds
            
            pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
            
        end
        
        for comp = 1:no_comparisons
        
            p_vals(comp, b, ch, 1) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
        
            p_vals(comp, b, ch, 2) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
            
        end
        
        p_vals = p_vals*all_dimensions(@sum, ~isnan(p_vals(:, b, :, :)));
        
        subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)pct_bp_high_for_test = max(pct_bp_high_for_test, eps);
        
        pct_bp_high_for_test = max(pct_bp_high_for_test, eps);
        
        if strcmp(measure, '_power')
            
            h = barwitherr(nanstd(pct_bp_high_for_test)*diag(1./sqrt(sum(~isnan(pct_bp_high_for_test)))), nanmean(pct_bp_high_for_test), 0.6);
            
            ylabel([chan_labels{ch}, 'Beta Power'])
            
            title([chan_labels{ch}, ', High ', band_labels{b}, ' Power (Per Trial)'])
            
        else
            
            h = barwitherr(nanstd(log(pct_bp_high_for_test))*diag(1./sqrt(sum(~isnan(log(pct_bp_high_for_test))))), nanmean(log(pct_bp_high_for_test)),...
                0.6, 'BaseValue', min(nanmean(log(pct_bp_high_for_test))) - 5);
            
            ylabel([chan_labels{ch}, 'Beta Density'])
            
            title([chan_labels{ch}, ', High ', band_labels{b}, ' Density (Per Trial)'])
            
        end
        
        box off
        
        bar_pos = get_bar_pos(h);
        
        for comp_type = 1:2
            
            bar_pairs = {};
            
            for comp = 1:no_comparisons
                
                if p_vals(comp, b, ch, comp_type) < .05
                    
                    bar_pairs = {bar_pairs{:}, [bar_pos(comparisons(comp, 1)), bar_pos(comparisons(comp, 2))]};
                    
                end
                
            end
            
            sigstar(bar_pairs, p_vals(p_vals(:, b, ch, comp_type) < .05, b, ch, comp_type))
            
        end
        
        xlim([0 (no_pds + 1)])
        
        set(gca, 'XTick', 1:no_pds, 'XTickLabel', pd_labels)
        
        % if b == 1
        
        % else
        % 
        %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
        % 
        % end
        
    end

    save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_', num2str(no_trials), 'trials_', short_band_labels{b}, '_barplot', measure])
    
end

%% Boxplots.
% 
% [p_greater, p_less] = deal(nan(no_bands, no_comparisons, no_chans));
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
%             
%         end
%         
%         for comp = 1:no_comparisons
%         
%             p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
%         
%             p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
%             
%         end
%         
%         subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         boxplot(pct_bp_high_for_test, 'labels', pd_labels)
%         
%         ylabel('High Beta Density')
%         
%         % if b == 1
%             
%         title([chan_labels{ch}, ', Boxplot of High ', band_labels{b}, ' Density (Per Trial)'])
%         
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials_', short_band_labels{b}, '_boxplot'])
%     
% end

%% Histograms.
% 
% [p_greater, p_less] = deal(nan(no_bands, no_comparisons, no_chans));
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
%             
%         end
%         
%         for comp = 1:no_comparisons
%         
%             p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
%         
%             p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
%             
%         end
%         
%         subplot(1, no_chans, ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         range = [all_dimensions(@min, pct_bp_high_for_test) all_dimensions(@max, pct_bp_high_for_test)];
%         
%         range = max(eps, range);
%         
%         [bin_edges, bin_centers] = make_bins(range(1), range(2), 50, 'log'); % 10, ''); %
%         
%         for pd = 1:no_pds
%             
%             [h, ~] = histc(pct_bp_high_for_test(:, pd), bin_edges);
%             
%             % h = max(eps, h);
%             
%             loglog(bin_centers, h(1:(end - 1))/sum(h(1:(end - 1))), pd_colors{pd})
%             
%             hold on
%             
%         end
%         
%         legend(pd_labels)
%         
%         xlabel('High Beta Density')
%         
%         ylabel('Proportion Observed')
%             
%         title([chan_labels{ch}, ', Histogram of High ', band_labels{b}, ' Density (Per Trial)'])
%         
%         % if b == 1
%         % 
%         %     title([chan_labels{ch}, ', ', band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials_', short_band_labels{b}, '_hist'])
%     
% end

end

function pos_bars = get_bar_pos(handle)

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end