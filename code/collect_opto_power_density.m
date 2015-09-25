function collect_opto_power_density(measure, norm_for_power, freqs, no_cycles, bands)

no_secs = 5; no_trials = 10;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
    BP_suffix = '';
    
else
    
    no_bands = size(bands, 1);
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

subject_matnames = {'st_m1_emxarch', 'st_m1_ali_post_carb_opto'};

channels = [1 2; 2 1];

no_folders = nan(1, 2);

for s = 1:2
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); % no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_pct_bp_high = nan(no_secs*no_trials, total_folders, no_pds, no_bands, 2);

for s = 1:2
   
    clear pct_bp_high
    
    load([subject_matnames{s}, BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials', measure, norm_for_power, '.mat'])
    
    if exist('BP_sec')
        
        pct_bp_high = BP_sec;
        
    end
    
    for ch = 1:2
    
        All_pct_bp_high(:, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :, ch) = pct_bp_high(:, :, :, :, channels(s, ch));
    
    end
    
end

clear pct_bp_high

pct_bp_high = All_pct_bp_high;

save(['OPTO', BP_suffix, '_pct_BP_high_laser_', num2str(no_trials), 'trials', norm_for_power, measure, '.mat'], 'pct_bp_high')