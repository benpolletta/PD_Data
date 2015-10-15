function collect_striatal_motor_power_density(measure, norm_for_power, freqs, no_cycles, bands)

epoch_secs = 150;

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
    BP_suffix = '';
    
else
    
    no_bands = size(bands, 1);
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

subject_matnames = {'st_m1', 'st_m1_ali'};

channels = [1 2; 2 1];

no_folders = nan(1, 2);

for s = 1:2
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_pct_bp_high = nan(epoch_secs, total_folders, no_pds, no_bands, no_chans);

for s = 1:2
   
    clear pct_bp_high
    
    load([subject_matnames{s}, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR', norm_for_power, measure, '.mat'])
    
    if exist('BP_sec', 'var')
        
        pct_bp_high = BP_sec;
        
    end
    
    for ch = 1:no_chans
    
        All_pct_bp_high(:, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :, ch) = pct_bp_high(:, :, :, :, channels(s, ch));
    
    end
    
end

clear pct_bp_high

pct_bp_high = All_pct_bp_high;

save(['STR_M1', BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', norm_for_power, measure, '.mat'], 'pct_bp_high')