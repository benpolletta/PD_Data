function PD_beta_blocks_rel_infusion_pre_post_PLV_plot(subject_mat, peak_suffix, epoch_secs, pd_handle, band_index_for_time, band_index_for_display, freqs, no_cycles, bands)

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
    
    PLV_name = [subj_mat_name, BP_suffix, '_pct_', short_band_labels{band_index_for_time}, '_high_',...
        num2str(epoch_secs/60), '_min_secs', pd_handle, '_PLV'];
    
else
    
    PLV_name = [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_',...
        short_band_labels{band_index_for_time}, pd_handle, '_PLV'];
    
    epoch_secs = str2double(pd_handle(2:(end - length('trials'))));
    
end

load([PLV_name, '.mat'])

no_freqs = sum(display_indices); 

no_dps = size(dP_sec, 2);

%% Group mean spectrum plots.

figure

array_names = {'Coh_sec', 'Coh_sec_pct', 'dP_sec', 'dP_sec_pct'};

long_array_names = {'Coherence', 'Coherence (%\Delta Baseline)', 'Phase of Coh.', 'Phase of Coh. (%\Delta Baseline)'};

no_arrays = length(array_names);

[rows, cols] = subplot_size(no_arrays);

figure

for a = 1:no_arrays
    
    eval(['PLV_data = ', array_names{a}, ';'])
    
    [PLV_mean, PLV_ci] = deal(nan(no_freqs, no_pds));
    
    for pd = 1:no_pds
        
        PLV_for_stats = reshape(PLV_data(display_indices, :, :, pd), no_freqs, no_folders*no_dps, 1);
        
        if strcmp(array_names{a}(1:2), 'dP_sec')
            
            PLV_mean(:, pd) = angle(nanmean(exp(sqrt(-1)*PLV_for_stats), 2));
    
            PLV_mean(PLV_mean(:, pd) < -pi/2, pd) = 2*pi + PLV_mean(PLV_mean(:, pd) < -pi/2, pd);
        
            PLV_ci(:, pd) = circ_confmean(PLV_for_stats, .05/bonferroni_count, [], [], 2);
            
        else
            
            PLV_mean(:, pd) = nanmean(PLV_for_stats, 2);
            
            PLV_ci(:, pd) = norminv(1 - .05/bonferroni_count, 0, 1)*nanstd(PLV_for_stats, [], 2)/sqrt(no_folders*no_dps);
            
        end
        
    end
    
    save([PLV_name, '_', array_names{a}, '_data_for_plot.mat'], 'PLV_mean', 'PLV_ci')
    
    subplot(rows, cols, a)
    
    if any(~isnan(PLV_mean))
        
        if any(~isnan(PLV_ci))
            
            boundedline(freqs(display_indices), PLV_mean, prep_for_boundedline(PLV_ci))
            
        else
            
            plot(freqs(display_indices), PLV_mean)
            
        end
        
    end
    
    axis tight
    
    title({[num2str(epoch_secs/60), ' Minutes, Densest ', band_labels{band_index_for_time}, ' Power']; long_array_names{a}})
    
    if a == 1
        
        legend(pd_labels)
        
    end
    
end

save_as_pdf(gcf, [PLV_name, '_', short_band_labels{band_index_for_display}])
