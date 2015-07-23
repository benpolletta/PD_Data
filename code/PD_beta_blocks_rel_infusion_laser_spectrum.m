function PD_beta_blocks_rel_infusion_laser_spectrum(subject_mat, norm, no_trials_analyzed, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3,7,length(freqs)); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = '';
    
else
    
    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2); 
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

pd_indices = cell(no_folders, 1);

trials = cell(no_folders, no_pds);

no_trials = nan(no_folders, no_pds);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    pd_indices{fo} = logical(laser_periods);
    
    for pd = 1:no_pds
        
        trials{fo, pd} = index_to_blocks(laser_periods(:, pd));
        
        no_trials(fo, pd) = size(trials, 1);
        
    end
    
end

if isempty(no_trials_analyzed)

    max_no_trials = all_dimensions(@max, no_trials);

else
    
    max_no_trials = no_trials_analyzed;
    
end
    
no_secs = max_no_trials*5;

WT_sec = nan(sum(band_indices{band_index}), no_secs, no_folders, no_pds, no_chans);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    [~, Spec_data] = get_BP(subj_name, outlier_lims(fo), norm, in_freqs, in_no_cycles, in_bands);
    
    if strcmp(norm, '')
        
        Spec_data = abs(Spec_data);
        
    end
    
    if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high_cum')
    
    else
        
        load([subj_name, BP_suffix, '_2sd_BP_high.mat'], 'BP_high')
        
        BP_high_cum = BP_high;
        
    end
    
    for ch = 1:no_chans
        
        for pd = 1:no_pds
            
            for tr = 1:min(max_no_trials, no_trials(fo, pd))
                
                trial_start = trials{fo, pd}(tr, 1);
                
                for sec = 1:5
                    
                    sec_start = trial_start + (sec - 1)*sampling_freq + 1;
                    
                    sec_end = trial_start + sec*sampling_freq;
                    
                    WT_sec(:, (tr - 1)*5 + sec, fo, pd, ch) = nanmean(Spec_data(sec_start:sec_end, band_indices{band_index}, ch))';
                    
                end
                
            end
            
        end
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_', short_band_labels{band_index},...
    '_', num2str(no_trials_analyzed), 'trials', norm, '_spectrum.mat'], 'WT_sec')

end
