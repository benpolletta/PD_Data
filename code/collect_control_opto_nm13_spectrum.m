function collect_control_opto_nm13_spectrum(peak_suffix, norm_for_power, band_index, freqs, no_cycles, bands)

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

subject_matnames = {'mice2_control_opto', 'st_m1_ali2_control_opto', 'st_m1_ali3_control_opto'};

no_mats = length(subject_matnames);

channels = repmat([2 1], 3, 1);

no_folders = nan(1, no_mats);

for s = 1:no_mats
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); % no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_WT_sec = nan(length(band_indices{band_index}), no_secs*no_trials, total_folders, no_pds, 2);

for s = 1:no_mats
   
    clear WT_sec
    
    load([subject_matnames{s}, BP_suffix, '_pct_', short_band_labels{band_index}, '_',...
        num2str(no_trials), 'trials', norm_for_power, '_spectrum.mat'])
    
    for ch = 1:2
        
        All_WT_sec(:, :, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, ch) = WT_sec(:, :, :, :, channels(s, ch));
        
    end
    
end

clear WT_sec

WT_sec = All_WT_sec;

save(['CONTROL_OPTO_no_mice3', BP_suffix, '_pct_', short_band_labels{band_index}, '_',...
    num2str(no_trials), 'trials', norm_for_power, '_spectrum.mat'], 'WT_sec')