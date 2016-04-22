function collect_control_opto_power_density(peak_suffix, measure, norm_for_power, freqs, no_cycles, bands)

no_secs = 5; no_trials = 10;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
    BP_suffix = peak_suffix;
    
else
    
    no_bands = size(bands, 1);
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

subject_matnames = {'mice2_control_opto', 'st_m1_ali2_control_opto', 'st_m1_ali3_control_opto'};

no_mats = length(subject_matnames);

channels = [2 1; 2 1; 2 1];

no_folders = nan(1, no_mats);

for s = 1:no_mats
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_pct_bp_high = nan(no_secs*no_trials, total_folders, no_pds, no_bands, no_chans);

for s = 1:no_mats
   
    clear pct_bp_high
    
    load([subject_matnames{s}, BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials', measure, norm_for_power, '.mat'])
    
    if exist('BP_sec')
        
        pct_bp_high = BP_sec;
        
    end
    
    for ch = 1:no_chans
    
        All_pct_bp_high(:, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :, ch) = pct_bp_high(:, :, :, :, channels(s, ch));
    
    end
    
end

clear pct_bp_high

pct_bp_high = All_pct_bp_high;

save(['CONTROL_OPTO_no_mice1_mice3', BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials', norm_for_power, measure, '.mat'], 'pct_bp_high')