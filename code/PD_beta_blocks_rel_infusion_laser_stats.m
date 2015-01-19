function PD_beta_blocks_rel_infusion_laser_stats(subject_mat)
    
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

no_secs = max_no_trials*5;

% high_type = {'', '_cum'}; no_types = length(high_type);

if isempty(dir([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials.mat']))

pct_bp_high = nan(no_secs, no_folders, no_pds, no_bands, no_chans);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    t = (1:size(BP_high_cum, 1))/sampling_freq;
    
    for ch = 1:no_chans
        
        for b = 1:no_bands
            
            figure(b)
            
            beta_blocks_plot = conv(BP_high_cum(:, b, ch), ones(sampling_freq, 1)/sampling_freq, 'same');
            
            subplot(no_folders, 2, (fo - 1)*2 + ch)
            
            plot(t/60, beta_blocks_plot)
            
            axis tight
            
            hold on
            
            handle = nan(no_pds, 1);
            
            for pd = 1:no_pds
                
                handle(pd) = plot(t(pd_indices{fo}(:, pd))/60, beta_blocks_plot(pd_indices{fo}(:, pd)), [pd_colors{pd}, '.']);
                
                for tr = 1:min(max_no_trials, no_trials(fo, pd))
                    
                    BP_trial = BP_high_cum(trials{fo, pd}(tr, 1):trials{fo, pd}(tr, 2), b, ch);
                    
                    BP_trial = reshape(BP_trial(1:5*sampling_freq), sampling_freq, 5);
                    
                    pct_bp_high((tr - 1)*5 + (1:5), fo, pd, b, ch) = nansum(BP_trial)'/sampling_freq;
                    
                end
                
            end
            
            axis tight
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
               
            elseif fo == no_folders
                
                xlabel('Time (Min.)')
                
            end
                    
            if ch == 1
            
                ylabel({folder; 'High Power Density per Sec.'})
                
                if fo == 1
                
                    legend(handle, pd_labels)
                
                end
                    
            end
            
        end
        
    end
    
end

for b = 1:no_bands
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials_', short_band_labels{b}])
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials.mat'], 'pct_bp_high')

else
    
load([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials.mat'])
    
end

no_comparisons = nchoosek(no_pds, 2);

comparisons = nchoosek(1:no_pds, 2);

[p_greater, p_less] = deal(nan(no_bands, no_comparisons, no_chans));

for b = 1:no_bands
    
    figure
    
    for ch = 1:no_chans
        
        pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
        
        for pd = 1:no_pds
            
            pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
            
        end
        
        for comp = 1:no_comparisons
        
            p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
        
            p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
            
        end
        
        subplot(1, no_chans, ch), % subplot(no_bands, 2, (b - 1)*2 + ch)
        
        boxplot(pct_bp_high_for_test, 'labels', pd_labels)
        
        ylabel('High Beta Density')
        
        % if b == 1
            
        title([chan_labels{ch}, ', Boxplot of High ', band_labels{b}, ' Density (Per Trial)'])
        
        % else
        % 
        %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
        % 
        % end
        
    end

    save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials_', short_band_labels{b}, '_boxplot'])
    
end

for b = 1:no_bands
    
    figure
    
    for ch = 1:no_chans
        
        pct_bp_high_for_test = nan(no_folders*no_secs, no_pds);
        
        for pd = 1:no_pds
            
            pct_bp_high_for_test(1:no_folders*no_secs, pd) = reshape(pct_bp_high(:, :, pd, b, ch), no_folders*no_secs, 1);
            
        end
        
        for comp = 1:no_comparisons
        
            p_greater(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'left');
        
            p_less(b, comp, ch) = ranksum(pct_bp_high_for_test(:, comparisons(comp, 1)), pct_bp_high_for_test(:, comparisons(comp, 2)), 'tail', 'right');
            
        end
        
        subplot(1, no_chans, ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
        
        range = [all_dimensions(@min, pct_bp_high_for_test) all_dimensions(@max, pct_bp_high_for_test)];
        
        range = max(eps, range);
        
        [bin_edges, bin_centers] = make_bins(range(1), range(2), 50, 'log'); % 10, ''); %
        
        for pd = 1:no_pds
            
            [h, ~] = histc(pct_bp_high_for_test(:, pd), bin_edges);
            
            % h = max(eps, h);
            
            loglog(bin_centers, h(1:(end - 1))/sum(h(1:(end - 1))), pd_colors{pd})
            
            hold on
            
        end
        
        legend(pd_labels)
        
        xlabel('High Beta Density')
        
        ylabel('Proportion Observed')
            
        title([chan_labels{ch}, ', Histogram of High ', band_labels{b}, ' Density (Per Trial)'])
        
        % if b == 1
        % 
        %     title([chan_labels{ch}, ', ', band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
        % 
        % else
        % 
        %     title([band_labels{b}, ', p-value = ', num2str(p_val(b, ch))])
        % 
        % end
        
    end

    save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_laser_trials_', short_band_labels{b}, '_hist'])
    
end

end