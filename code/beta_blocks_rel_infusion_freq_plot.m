function beta_blocks_rel_infusion_freq_plot(subject_mat, norm, hist_norm, band_index, epoch_secs, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

if ~isempty(epoch_secs)
    
    epoch_label = sprintf('_%.1f_min', epoch_secs/60);
    
else
    
    epoch_label = '';
    
end

close('all')

subject_mat_name = subject_mat(1:(end - length('_subjects.mat')));

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, BP_suffix, '_wt.mat'], 'sampling_freq')

cycle_lengths = sampling_freq./freqs;

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

freq_lims = [min(freqs(band_indices{band_index})) max(freqs(band_indices{band_index}))];

no_pds = length(pd_labels);

no_chans = length(chan_labels);

%% Group average plot (by datapoint). 

[Freq_hist, Spec_high_beta_mean, Spec_high_beta_std] = deal(nan(sum(band_indices{band_index}), no_pds, no_chans));
    
for ch = 1:no_chans
    
    Freq_high_beta = cell(no_pds, 1);
    
    for pd = 1:no_pds
        
        Freq_data = load([subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index},...
            '_freqs', norm, '_ch', num2str(ch), '_', pd_labels{pd}, '.txt']);
        
        Freq_high_beta{pd} = Freq_data(:, 1);
        
        [h, ~] = hist(Freq_high_beta{pd}, freqs(band_indices{band_index}));
        
        if strcmp(hist_norm, '_cycles')
        
            Freq_hist(:, pd, ch) = h./cycle_lengths(band_indices{band_index});
        
        else
        
            Freq_hist(:, pd, ch) = h/sum(h);
        
        end
        
        figure(1)
        
        subplot(1, no_chans, ch)
        
        plot(freqs(band_indices{band_index}), h/sum(h), pd_colors{pd}, 'LineWidth', 2)
        
        hold on
        
        Spec_high_beta_mean(:, pd, ch) = nanmean(Freq_data(:, 2:end))';
        
        Spec_high_beta_std(:, pd, ch) = nanstd(Freq_data(:, 2:end))'/sqrt(size(Spec_high_beta_std, 1)/(sampling_freq/8));
        
    end
    
    axis tight
    
    xlabel('Freq. (Hz)')
    
    ylabel('Proportion High Beta Datapoints Observed')
    
    if ch == 1
    
        legend(pd_labels)
    
    end
    
    title([chan_labels{ch}, ', Histogram of High Beta Frequencies'])
    
    if no_pds >= 2
        
        comparisons = nchoosek(1:no_pds, 2);
        
        comparisons(~any(comparisons == 2, 2), :) = [];
        
        no_comps = size(comparisons, 1);
        
        comp_colors = colormap(hsv(no_comps + 1));
        
        [handles, p_vals] = deal(nan(no_comps, 1));
        
        comp_labels = cell(no_comps, 1);
        
        % comp_titles = cell(no_comps, 1);
        
        figure(2)
        
        subplot(1, no_chans, ch)
        
        for comp = 1:no_comps
            
            h = qqplot(Freq_high_beta{comparisons(comp, 1)}, Freq_high_beta{comparisons(comp, 2)});
            
            hold on
            
            plot(freq_lims, freq_lims, 'k')
                
            [~, p] = kstest2(Freq_high_beta{comparisons(comp, 1)}, Freq_high_beta{comparisons(comp, 2)});
            
            p_vals(comp) = p;
            
            if no_pds == 2
            
                xlabel(pd_labels{1}), ylabel(pd_labels{2})
        
                title({[chan_labels{ch}, ' Pre- vs. Post-Infusion Frequency Distribution']; ['p = ', num2str(p), ', Kolmogorov-Smirnov Test']})
            
            else
                
                handles(comp) = h(1);
                
                comp_labels{comp} = [pd_labels{comparisons(comp, 1)}, ' vs. ', pd_labels{comparisons(comp, 2)}];
                
                comp_labels{comp} = [comp_labels{comp}, ', p = ', num2str(p*no_comps), ', K-S Test'];
                
                set(h(1), 'Marker', 'x', 'MarkerEdgeColor', comp_colors(comp, :))
                
                set(h(2), 'Color', comp_colors(comp, :))
                
                set(h(3), 'Color', comp_colors(comp, :))
                
            end
            
        end
            
        xlim(freq_lims), ylim(freq_lims)
        
        if no_comps > 2
            
            legend(handles, comp_labels)
            
            title(chan_labels{ch}) % title(comp_titles)
            
        end
            
    end
    
end

save([subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index}, '_freqs',...
    norm, hist_norm, '_group_stats'], 'Freq_hist', 'Spec_high_beta_mean', 'Spec_high_beta_std')

save_as_pdf(1, [subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index},...
    '_freqs', norm, hist_norm, '_group_hist'])

if no_pds >= 2
    
    save_as_pdf(2, [subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index},...
        '_freqs', norm, hist_norm, '_group_qq'])
    
end

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)

    boundedline(freqs(band_indices{band_index}), Spec_high_beta_mean(:, :, ch), prep_for_boundedline(Spec_high_beta_std(:, :, ch)), 'cmap', pd_cmap)
      
    axis tight

    xlabel('Freq. (Hz)')
    
    ylabel('Mean \pm S.D. Power')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Mean Power During High Beta'])
    
end

save_as_pdf(gcf, [subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index}, '_freqs', norm, hist_norm, '_group_mean'])

%% Group average plot (by recording).

[Freq_hist, Spec_high_beta_mean, Spec_high_beta_std] = deal(nan(no_folders, sum(band_indices{band_index}), no_pds, no_chans));

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index}, '_freqs', norm, '.mat'])
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
        
        load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
        
        pd_indices = laser_periods;
        
    else
        
        pd_indices = nan(length(t), 2);
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    end
    
    pd_indices = logical(pd_indices);
    
    for ch = 1:no_chans
        
        for pd = 1:no_pds
            
            [h, ~] = hist(Freqs_high_beta(pd_indices(:, pd), ch), freqs(band_indices{band_index}));
            
            Freq_hist(fo, :, pd, ch) = h'/sum(h);
            
            Spec_high_beta_mean(fo, :, pd, ch) = nanmean(Spec_high_beta(pd_indices(:, pd), :, ch));
            
            Spec_high_beta_std(fo, :, pd, ch) = nanstd(Spec_high_beta(pd_indices(:, pd), :, ch));
            
        end
        
    end
    
end

save([subject_mat_name, BP_suffix, epoch_label, '_beta_block_', short_band_labels{band_index}, '_freqs',...
    norm, hist_norm, '_individual_stats'], 'Freq_hist', 'Spec_high_beta_mean', 'Spec_high_beta_std')

Freq_hist_mean = permute(nanmean(Freq_hist), [2 3 4 1]);

Freq_hist_std = permute(nanstd(Freq_hist), [2 3 4 1])/sqrt(no_folders);

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)
    
    boundedline(freqs(band_indices{band_index})', Freq_hist_mean(:, :, ch),...
        prep_for_boundedline(Freq_hist_std(:, :, ch)), 'cmap', pd_cmap)

    axis tight
    
    xlabel('Freq. (Hz)')
    
    xlim(freq_lims)
    
    ylabel('Proportion of High Beta')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Histogram of High Beta Frequency'])
   
end

save_as_pdf(gcf, [subject_mat_name, BP_suffix, epoch_label, '_beta_block_freqs', norm, '_individual_hist'])

Spec_high_beta_mean_plot = permute(nanmean(Spec_high_beta_mean), [2 3 4 1]);

Spec_high_beta_std_plot = permute(nanstd(Spec_high_beta_mean), [2 3 4 1])/sqrt(no_folders);

figure

for ch = 1:no_chans
    
    subplot(1, no_chans, ch)
    
    boundedline(freqs(band_indices{band_index})', Spec_high_beta_mean_plot(:, :, ch),...
        prep_for_boundedline(Spec_high_beta_std_plot(:, :, ch)), 'cmap', pd_cmap)
    
    xlim(freq_lims)

    axis tight
    
    xlabel('Freq. (Hz)')
    
    ylabel('Mean \pm S.D. Power')
    
    if ch == 1
        
        legend(pd_labels)
        
    end
    
    title([chan_labels{ch}, ', Mean Power During High Beta'])
   
end

save_as_pdf(gcf, [subject_mat_name, BP_suffix, epoch_label, '_beta_block_freqs', norm, hist_norm, '_individual_mean'])