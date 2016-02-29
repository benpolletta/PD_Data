function PD_beta_blocks_rel_infusion_pre_post_spectrum_plot(subject_mat, peak_suffix, epoch_secs, pd_handle, norm, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)

% Leave epoch_secs empty when using for optogenetics data, and enter
% '_ntrials' for the argument pd_handle.

load('bonferroni_count.mat')

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

display_indices = band_indices{band_index_for_display};

no_chans = length(chan_labels);

no_pds = length(pd_labels);

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));

if ~isempty(epoch_secs)
    
    spectrum_name = [subj_mat_name, BP_suffix, '_pct_', short_band_labels{band_index_for_time}, '_high_',...
        num2str(epoch_secs/60), '_min_secs', pd_handle, norm, '_spectrum'];
    
else
    
    spectrum_name = [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_',...
        short_band_labels{band_index_for_time}, pd_handle, norm, '_spectrum'];
    
    epoch_secs = str2double(pd_handle(2:(end - length('trials'))));
    
end

load([spectrum_name, '.mat'])

no_freqs = sum(display_indices); 

no_dps = size(WT_sec, 2);

%% Group mean spectrum plots.

figure
        
% matlabpool open; opt = statset('UseParallel', true);

for ch = 1:no_chans
    
    [WT_mean, WT_ci] = deal(nan(no_freqs, no_pds)); % [WT_mean, WT_se] [WT_median, WT_median_low, WT_median_high] 
    
    for pd = 1:no_pds
        
        WT_for_stats = reshape(WT_sec(display_indices, :, :, pd, ch), no_freqs, no_folders*no_dps, 1);
        
        WT_mean(:, pd) = nanmean(WT_for_stats, 2); % .*freq_multiplier, 2);
        
        WT_ci(:, pd) = norminv(1 - .05/bonferroni_count, 0, 1)*nanstd(WT_for_stats, [], 2)/sqrt(no_folders*no_dps); % .*freq_multiplier, [], 2)/sqrt(epoch_secs);
        
        % WT_median(:, pd) = nanmedian(WT_for_stats, 2);
        % 
        % WT_median_distn = jackknife(@nanmedian, WT_for_stats', 'Options', opt);
        % 
        % WT_median_low(:, pd) = 1.5*range(WT_median_distn)'; % norminv(1 - .25/(no_freqs*7*2*3), 0, 1)*std(WT_median_distn)'; % quantile(WT_median_distn, .025)';
        % 
        % WT_median_high(:, pd) = 1.5*range(WT_median_distn)'; % norminv(.25/(no_freqs*7*2*3), 0, 1)*std(WT_median_distn)'; % quantile(WT_median_distn, .975)';
        
        % WT_median_low(:, pd) = quantile(WT_for_mean, .25, 2);
        % 
        % WT_median_high(:, pd) = quantile(WT_for_mean, .75, 2);
        
    end
    
    save([spectrum_name, '_ch', num2str(ch), '_data_for_plot.mat'], 'WT_mean', 'WT_ci')
    
    % WT_median_low = WT_median - WT_median_low;
    % 
    % WT_median_high = WT_median_high - WT_median;
    
    subplot(1, no_chans, ch)
    
    if strcmp(norm, '_pct')
        
        boundedline(freqs(display_indices), WT_mean, prep_for_boundedline(WT_ci))
        
        % boundedline(freqs(display_indices), WT_median, prep_for_boundedline(WT_median_low, WT_median_high))
        
    else
        
        freq_multiplier = repmat(freqs(display_indices)', 1, no_pds);
        
        boundedline(freqs(display_indices), WT_mean.*freq_multiplier, prep_for_boundedline(WT_ci.*freq_multiplier))
        
    end
    
    axis tight
    
    title({[chan_labels{ch}, ', ', num2str(epoch_secs/60), ' Minutes,'];['Densest High Power, ', band_labels{band_index_for_time}]})
    
    legend(pd_labels)
    
end

% matlabpool close

save_as_pdf(gcf, [spectrum_name, '_', short_band_labels{band_index_for_display}])
