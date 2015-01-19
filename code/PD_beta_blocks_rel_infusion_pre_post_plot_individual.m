function PD_beta_blocks_rel_infusion_pre_post_plot_individual(subject_mat, epoch_secs, test_handle, pd_handle)
    
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

% pd_colors = {'g', 'r'};

% norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);
% 
% long_norms = {'', ', Increase Over Baseline Power', ', % Total Power', ', Increase in % Total Power Over Baseline'};

% high_type = {'', '_cum'}; no_types = length(high_type);

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));
    
load([subj_mat_name, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

if exist('BP_sec')
    
    pct_bp_high = BP_sec;
    
end

%% Individual bar plots.

for b = 1:no_bands
    
    figure
    
    for fo = 1:no_folders
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        p_val = nan(no_bands, no_chans);
        
        for ch = 1:no_chans
            
            pct_bp_high_for_test = nan(epoch_secs, no_pds);
            
            for pd = 1:no_pds
                
                pct_bp_high_for_test(:, pd) = pct_bp_high(:, fo, pd, b, ch);
                
            end
            
            if no_pds == 2 %for comp = 1:no_comparisons
                
                if strcmp(test_handle, 't-test')
                
                    [~, p_val(b, ch)] = ttest(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 2), [], 'left'); % 'tail', 'left');
                
                elseif strcmp(test_handle, 'ranksum')
                    
                    p_val(b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 2), 'tail', 'left');
                    
                end
                
            end
            
            % p_val = p_val*all_dimensions(@sum, ~isnan(p_val))*no_folders;
            
            subplot(no_chans, no_folders, (ch - 1)*no_folders + fo), % subplot(no_bands, 2, (b - 1)*2 + ch)
            
            pct_bp_high_for_test = max(pct_bp_high_for_test, eps);
            
            if exist('BP_sec')
                
                h = barwitherr(nanstd(pct_bp_high_for_test)*diag(1./sqrt(sum(~isnan(pct_bp_high_for_test)))), nanmean(pct_bp_high_for_test), 0.6);
                
                if fo == 1
                    
                    ylabel({chan_labels{ch}; 'Beta Power'})
                    
                end
                
            else
            
                h = barwitherr(nanstd(log(pct_bp_high_for_test))*diag(1./sqrt(sum(~isnan(log(pct_bp_high_for_test))))), nanmean(log(pct_bp_high_for_test)),...
                    0.6, 'BaseValue', min(nanmean(log(pct_bp_high_for_test))) - 5);
                
                if fo == 1
                    
                    ylabel({chan_labels{ch}; 'Beta Density'})
                    
                end
            
            end
                
            box off
            
            if p_val(b, ch) < 0.05
                
                bar_pos = get_bar_pos(h);
                
                sigstar(bar_pos, p_val(b, ch))
                
            end
            
            if fo == 1
            
                title({[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
                    [folder, ', ', band_labels{b}];['p-value = ', num2str(p_val(b, ch)), ' (t-test)']})
                
            else
                
                title({[folder, ','];['p-value = ', num2str(p_val(b, ch)), ' (', test_handle, ')']})
            
            end
            
            xlim([0 3])
            
            set(gca, 'XTickLabel', pd_labels)
            
        end
    
        save([subj_name, '_pct_BP_high_pvals'], 'p_val')
        
    end
        
    save_as_pdf(gcf, [subj_mat_name, '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_barplot_individual_', test_handle, pd_handle])
    
end

%% Group boxplots.
%
% p_val = nan(no_bands, no_chans);
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*epoch_secs, no_pds);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(:, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*epoch_secs, 1);
%             
%         end
%         
%         if no_pds == 2 %for comp = 1:no_comparisons
%         
%             p_val(b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 2), 'tail', 'left');
%             
%         end
%         
%         subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         boxplot(pct_bp_high_for_test, 'labels', pd_labels)
%         
%         ylabel('High Beta Density')
%         
%         % if b ==
%             
%             title({[chan_labels{ch}, ','];[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
%                 [band_labels{b}, ', p-value = ', num2str(p_val(b, ch)), ' (ranksum test)']})
%         
%         % else
%         % 
%         %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
%         % 
%         % end
%         
%     end
% 
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_boxplot'])
%     
% end

%% Group histograms.
% 
% p_val = nan(no_bands, no_chans);
% 
% for b = 1:no_bands
%     
%     figure
%     
%     for ch = 1:no_chans
%         
%         pct_bp_high_for_test = nan(no_folders*epoch_secs, 2);
%         
%         for pd = 1:no_pds
%             
%             pct_bp_high_for_test(:, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*epoch_secs, 1);
%             
%         end
%         
%         if no_pds == 2 % for comp = 1:no_comparisons
%         
%             p_val(b, ch) = ranksum(pct_bp_high_for_test(:, 1), pct_bp_high_for_test(:, 1), 'tail', 'left');
%             
%         end
%         
%         subplot(1, no_chans, ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
%         
%         range = [all_dimensions(@min, pct_bp_high_for_test) all_dimensions(@max, pct_bp_high_for_test)];
%         
%         range = max(eps, range);
%         
%         [bin_edges, bin_centers] = make_bins(range(1), range(2), 100, 'log'); % 25, ''); % 
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
%         title({[chan_labels{ch}, ','];[num2str(epoch_secs/60), ' Minutes of Densest High Power'];...
%             [band_labels{b}, ', p-value = ', num2str(p_val(b, ch)), ' (ranksum test)']})
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
%     save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{b}, '_hist'])
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