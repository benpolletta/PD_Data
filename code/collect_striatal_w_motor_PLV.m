function collect_striatal_w_motor_PLV(peak_suffix, epoch_secs, band_index, freqs, no_cycles, bands)

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

[All_Coh_sec, All_Coh_pct_sec, All_dP_sec, All_dP_pct_sec] = deal(nan(no_freqs, epoch_secs, total_folders, no_pds));

array_names = {'Coh', 'Coh_pct', 'dP', 'dP_pct'};

for s = 1:no_groups
   
    clear WT_sec
    
    load([subject_matnames{s}, BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
        num2str(epoch_secs/60), '_min_secs_by_STR_PLV.mat'])
    
    All_dP_sec(:, :, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :) = ((-1)^channels(s))*dP_sec;
    
    for a = [1 2 4]
        
        eval(['All_', array_names{a}, '_sec(:, :, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1)) = ', array_names{a}, '_sec;'])
        
    end
    
end

% clear Coh_sec Coh_pct_sec dP_sec dP_pct_sec

for a = 1:length(array_names)
    
    clear(array_names{a})
   
    eval([array_names{a}, '_sec = All_', array_names{a}, '_sec;'])
    
end

save(['STR_M1', BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
    num2str(epoch_secs/60), '_min_secs_PLV.mat'], 'Coh_sec', 'Coh_pct_sec', 'dP_sec', 'dP_pct_sec')