function PD_beta_blocks_rel_infusion_pre_post_spectrum(subject_mat, epoch_secs, pd_handle, norm, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders); no_pds = length(pd_labels); no_chans = length(chan_labels);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

WT_sec = nan(sum(band_indices{band_index}), epoch_secs, no_folders, no_pds, no_chans);
    
load([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    [~, Spec_data] = get_BP(subj_name, outlier_lims(fo), norm, freqs, no_cycles, bands);
    
    if strcmp(norm, '')
        
        Spec_data = abs(Spec_data);
        
    end
    
    t = (1:size(Spec_data, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    clear pd_indices
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), 1);
        
    end
    
    pd_indices = logical(pd_indices);
    
    for ch = 1:no_chans
        
        for pd = 1:no_pds
            
            bp_max_start = All_bp_max_start(fo, ch, band_index, pd);
           
            for sec = 1:epoch_secs
                
                sec_start = max(bp_max_start + (sec - 1)*sampling_freq + 1, 1);
                
                sec_end = min(max(bp_max_start + sec*sampling_freq, 1), find(pd_indices(:, pd) == 1, 1, 'last'));
                
                WT_sec(:, sec, fo, pd, ch) = nanmean(Spec_data(sec_start:sec_end, band_indices{band_index}, ch))';
                
            end
            
        end
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_', short_band_labels{band_index}, '_high_',...
    num2str(epoch_secs/60), '_min_secs', pd_handle, norm, '_spectrum.mat'], 'freqs', 'band_indices', 'band_index', 'WT_sec')

end