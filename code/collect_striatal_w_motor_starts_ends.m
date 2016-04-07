function collect_striatal_w_motor_starts_ends(peak_suffix, freqs, no_cycles, bands)

epoch_secs = 150;

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

subject_matnames = {'st_m1', 'st_m1_ali', 'st_m1_ali2'};

channels = [1 2 2];

no_groups = length(subject_matnames);

no_folders = nan(1, no_groups);

for s = 1:no_groups
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); % no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

[bp_max_start, bp_max_end] = deal(nan(total_folders, 1, no_bands, no_pds));

for s = 1:no_groups
   
    clear All_bp_max_start All_bp_max_end
    
    load([subject_matnames{s}, BP_suffix, '_pct_BP_high_',...
        num2str(epoch_secs/60), '_min_secs_by_STR.mat'])
    
    bp_max_start((subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :, :) =...
        All_bp_max_start(:, channels(s), :, :);
    
    bp_max_end((subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :, :, :) =...
        All_bp_max_end(:, channels(s), :, :);
    
end

clear All_bp_max_start All_bp_max_end

All_bp_max_start = bp_max_start; All_bp_max_end = bp_max_end;

save(['STR_M1', BP_suffix, '_pct_BP_high_',...
    num2str(epoch_secs/60), '_min_secs_by_STR.mat'], 'All_bp_max_start', 'All_bp_max_end')