function collect_striatal_spectrum(epoch_secs, norm_for_power, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end
    
no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

subject_matnames = {'st_m1', 'st_stn'};

channels = [1 2];

no_folders = nan(1, 2);

for s = 1:2
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
end

no_pds = length(pd_labels); % no_chans = length(chan_labels);

subj_mat_limits = [0 cumsum(no_folders)];

total_folders = subj_mat_limits(end);

All_WT_sec = nan(sum(band_indices{band_index}), epoch_secs, total_folders, no_pds);

for s = 1:2
   
    clear WT_sec
    
    load([subject_matnames{s}, BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
        num2str(epoch_secs/60), '_min_secs_by_STR', norm_for_power, '_spectrum.mat'])
    
    All_WT_sec(:, :, (subj_mat_limits(s) + 1):subj_mat_limits(s + 1), :) = WT_sec(:, :, :, :, channels(s));
    
end

clear WT_sec

WT_sec = All_WT_sec;

save(['STR', BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
    num2str(epoch_secs/60), '_min_secs', norm_for_power, '_spectrum.mat'], 'WT_sec')