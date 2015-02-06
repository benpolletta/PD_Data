function combine_carb_6OHDA_power_density(epoch_secs, measure, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    no_bands = 6;
    
    BP_suffix = '';
    
else
    
    no_bands = size(bands, 1);
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
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

All_pct_bp_high = nan(epoch_secs, max_folders, 3, no_bands, 2);

for s = 1:no_subj_mats
   
    clear pct_bp_high
    
    load([subject_matnames{s}, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handles{s}, measure, '.mat'])
    
    if exist('BP_sec', 'var')
        
        pct_bp_high = BP_sec;
        
    end
    
    no_folders = size(pct_bp_high, 2);
    
    for fo = 1:no_folders
        
        All_pct_bp_high(:, fo, in_periods{s}, :, in_channels{s}(fo, :)) = pct_bp_high(:, fo, out_periods{s}, :, out_channels{s}(fo, :));
        
    end
    
end

clear pct_bp_high

pct_bp_high = All_pct_bp_high;

save(['CARB_6OHDA', BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', measure, '.mat'], 'pct_bp_high')