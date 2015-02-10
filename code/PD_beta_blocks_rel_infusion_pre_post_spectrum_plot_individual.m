function PD_beta_blocks_rel_infusion_pre_post_spectrum_plot_individual(subject_mat, epoch_secs, pd_handle, norm, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

% if exist([folders{1}, '/', prefixes{1}, '_wt.mat'])
%     
%     load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')
%     
% else
%     
%     sampling_freq = 500;
%     
% end

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
   
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
    
load([subj_mat_name, BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
    num2str(epoch_secs/60), '_min_secs', pd_handle, norm, '_spectrum.mat'])

no_freqs = size(WT_sec, 1);

% no_comparisons = nchoosek(no_pds, 2);
%     
% p_vals = nan(size(WT_sec, 1), no_comparisons, no_folders, no_bands, no_chans, 2);
%     
% All_p_vals = nan(size(WT_sec, 1), no_comparisons, no_bands, no_chans, 2);

%% Individual bar plots.

figure

All_mean = nan(no_freqs, no_folders, no_pds, no_chans);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    for ch = 1:no_chans
        
        [WT_mean, WT_se] = deal(nan(no_freqs, no_pds));
        
        for pd = 1:no_pds
            
            WT_mean(:, pd) = nanmean(WT_sec(:, :, fo, pd, ch), 2);
            
            WT_se(:, pd) = nanstd(WT_sec(:, :, fo, pd, ch), [], 2)/sqrt(epoch_secs);
            
        end
        
        % subj_p_vals = run_stats(pct_bp_high_for_test, test_handle);
        
        % p_vals(:, fo, b, ch, :) = permute(subj_p_vals, [1 3 4 5 2]);
        
        % p_val = p_val*all_dimensions(@sum, ~isnan(p_val))*no_folders;
        
        subplot(no_chans, no_folders, (ch - 1)*no_folders + fo), % subplot(no_bands, 2, (b - 1)*2 + ch)
        
        boundedline(freqs(band_indices{band_index}), log(WT_mean), prep_for_boundedline(log(WT_se)))
        
        axis tight
        
        % plot_data(pct_bp_high_for_test, fo, exist('BP_sec', 'var') || strcmp(pd_handle, '_power'), subj_p_vals)
        
        All_mean(:, fo, :, ch) = permute(WT_mean, [1 3 2]);
        
        if fo == 1
            
            title({[num2str(epoch_secs/60), ' Minutes of Densest High Power'];[folder, ', ', band_labels{b}]})
            
            legend(pd_labels)
            
        else
            
            title(folder)
            
        end
        
    end
    
end

save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{band_index},...
    '_spectrum_individual', pd_handle, norm])

%% Stats & figures treating each individual as an observation.

figure

for ch = 1:no_chans
    
    % across_p_vals = run_stats(All_mean(:, :, ch), test_handle);
    % 
    % All_p_vals(:, b, ch, :) = permute(across_p_vals, [1 3 4 2]);
    
    for pd = 1:no_pds
        
        All_mean_mean(:, pd) = nanmean(All_mean(:, :, pd, ch), 2);
        
        All_mean_se(:, pd) = nanstd(All_mean(:, :, pd, ch), [], 2)/sqrt(no_folders);
        
    end
    
    subplot(1, no_chans, ch)
    
    boundedline(freqs(band_indices{band_index}), log(All_mean_mean), prep_for_boundedline(log(All_mean_se)))
    
    axis tight
    
    title([chan_labels{ch}, ', ', num2str(epoch_secs/60), ' Minutes of Densest High Power, ', band_labels{band_index}])
    
    legend(pd_labels)
    
end

save_as_pdf(gcf, [subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_mins_', short_band_labels{band_index},...
    '_spectrum_individual_avg', pd_handle, norm])

end