function collect_striatal_freq(epoch_secs, measure)

subject_matnames = {'st_m1', 'st_stn'};

channels = [1 2];

no_folders = nan(1, 2);

for s = 1:2
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_pct_bp_high = nan(epoch_secs, total_folders, no_pds, 6);

for s = 1:2
   
    clear pct_bp_high
    
    load([subject_matnames{s}, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR', measure, '.mat'])
    
    if exist('BP_sec')
        
        pct_bp_high = BP_sec;
        
    end
    
    All_pct_bp_high(:, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :) = pct_bp_high(:, :, :, :, channels(s));
    
end

clear pct_bp_high

pct_bp_high = All_pct_bp_high;

save(['STR_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', measure, '.mat'], 'pct_bp_high')