function collect_opto_spectrum(peak_suffix, norm_for_power, band_index, freqs, no_cycles, bands)

no_secs = 5; no_trials = 10;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = peak_suffix;
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

no_freqs = length(freqs);
    
no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

subject_matnames = {'st_m1_emxarch', 'st_m1_ali_post_carb_opto', 'st_m1_ali2_post_carb_opto'};

no_mats = length(subject_matnames);

channels = [1 2; 2 1; 2 1];

no_folders = nan(1, no_mats);

for s = 1:no_mats
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); % no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

array_names = {'WT_sec', 'WT_trial_normed_sec', 'WT_trial_normed_mean', 'WT_trial_normed_se'};

for a = 1:2

    eval(sprintf('All_%s = nan(length(freqs), no_secs*no_trials, total_folders, no_pds, 2);', array_names{a}))

end

for a = 3:4
   
    eval(sprintf('All_%s = nan(length(freqs), total_folders, no_pds, 2);', array_names{a}))
    
end

for s = 1:no_mats
   
    clear WT_sec
    
    load([subject_matnames{s}, BP_suffix, '_pct_', short_band_labels{band_index}, '_',...
        num2str(no_trials), 'trials', norm_for_power, '_spectrum.mat'])
    
    for ch = 1:2
        
        for a = 1:2
        
            eval(sprintf('All_%s(:, :, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, ch) = %s(:, :, :, :, channels(s, ch));', array_names{a}, array_names{a}))
        
        end
        
        for a = 3:4
        
            eval(sprintf('All_%s(:, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, ch) = %s(:, :, :, channels(s, ch));', array_names{a}, array_names{a}))
        
        end
        
    end
    
end

clear WT_sec WT_trial_normed_sec WT_trial_normed_mean WT_trial_normed_se

for a = 1:length(array_names)

    eval(sprintf('%s = All_%s;', array_names{a}, array_names{a}))

end

save(['OPTO', BP_suffix, '_pct_', short_band_labels{band_index}, '_',...
    num2str(no_trials), 'trials', norm_for_power, '_spectrum.mat'], 'WT_sec', 'WT_trial_normed_sec', 'WT_trial_normed_mean', 'WT_trial_normed_se')