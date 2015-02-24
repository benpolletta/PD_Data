function combine_carb_6OHDA_spectrum(epoch_secs, norm_for_power, band_index, freqs, no_cycles, bands)

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

subject_matnames = {'STR', 'st_m1', 'st_m1_6OHDA'};

no_subj_mats = length(subject_matnames);

pd_handles = {'', '_by_STR', ''};

no_folders = nan(1, 2);

[in_channels, out_channels, in_periods, out_periods] = deal(cell(no_subj_mats, 1));

for s = 1:no_subj_mats
    
    load([subject_matnames{s}, '_subjects.mat'])
    
    no_folders(s) = length(folders);
    
    if exist('chan_id', 'var')
        
        out_channels{s} = [chan_id; 3-chan_id]';
        
        in_channels{s} = cumsum(ones(size(out_channels{s})), 2);
        
        in_periods{s} = 3;
        
        out_periods{s} = 1;
        
    else
        
        out_channels{s} = s*ones(length(folders), 1);
    
        in_channels{s} = s*ones(length(folders), 1);
        
        in_periods{s} = [1 2];
        
        out_periods{s} = [1 2];
        
    end
        
end

max_folders = max(no_folders); clear no_folders

All_WT_sec = nan(sum(band_indices{band_index}), epoch_secs, max_folders, 3, 2);

for s = 1:no_subj_mats
   
    clear WT_sec
    
    load([subject_matnames{s}, BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
        num2str(epoch_secs/60), '_min_secs', pd_handles{s}, norm_for_power, '_spectrum.mat'])
    
    no_folders = size(WT_sec, 3);
    
    for fo = 1:no_folders
        
        All_WT_sec(:, :, fo, in_periods{s}, in_channels{s}(fo, :)) = WT_sec(:, :, fo, out_periods{s}, out_channels{s}(fo, :));
        
    end
    
end

clear WT_sec

WT_sec = All_WT_sec;

save(['CARB_6OHDA', BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
    num2str(epoch_secs/60), '_min_secs', norm_for_power, '_spectrum.mat'], 'WT_sec')