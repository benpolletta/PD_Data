function PD_beta_blocks_rel_infusion_laser_PLV(subject_mat, peak_suffix, no_trials_analyzed, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3,7,length(freqs)); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = peak_suffix;
    
else
    
    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, BP_suffix, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2); 
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

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
        
        no_trials(fo, pd) = size(trials{fo, pd}, 1);
        
    end
    
end

if isempty(no_trials_analyzed)

    max_no_trials = all_dimensions(@max, no_trials);

else
    
    max_no_trials = no_trials_analyzed;
    
end
    
no_secs = max_no_trials*5;

[Coh_sec, Coh_sec_pct, dP_sec, dP_sec_pct] = deal(nan(sum(band_indices{no_bands}), no_secs, no_folders, no_pds));

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    folder_trial_no = min(max_no_trials, no_trials(fo, pd));
    
    subj_name = [folder,'/',prefix];
    
    [~, Spec_data] = get_BP(subj_name, peak_suffix, outlier_lims(fo), '', in_freqs, in_no_cycles, in_bands);
    
    Phase_data = angle(Spec_data);
    
    dPhase_data = diff(Phase_data, [], 3);
    
    for pd = 1:no_pds
        
        dP_trials = nan(folder_trial_no*5*sampling_freq, sum(band_indices{no_bands}));
        
        for tr = 1:folder_trial_no
            
            trial_start = trials{fo, pd}(tr, 1);
            
            trial_end = trial_start + 5*sampling_freq - 1;
            
            dP_trials((tr - 1)*5*sampling_freq + (1:5*sampling_freq), :) = dPhase_data(trial_start:trial_end, band_indices{no_bands});
            
        end
        
        dP_trials = nans_to_end(dP_trials);
        
        for sec = 1:folder_trial_no*5
            
            sec_start = (sec - 1)*sampling_freq + 1;
            
            sec_end = sec*sampling_freq;
            
            MRV_sec = nanmean(exp(sqrt(-1)*dP_trials(sec_start:sec_end, :)))';
            
            Coh_sec(:, sec, fo, pd) = abs(MRV_sec);
            
            dP_sec(:, sec, fo, pd) = angle(MRV_sec);
            
        end
        
        Coh_baseline = nanmean(Coh_sec(:, :, fo, pd), 2);
        
        dP_baseline = angle(nanmean(exp(sqrt(-1)*dP_sec(:, :, fo, pd)), 2));
        
        Coh_sec_pct(:, :, fo, pd) = (100*(Coh_sec(:, :, fo, pd)')*diag(1./Coh_baseline) - 100)';
        
        dP_sec_pct(:, :, fo, pd) = (100*(dP_sec(:, :, fo, pd)')*diag(1./dP_baseline) - 100)';
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_', short_band_labels{no_bands},...
    '_', num2str(no_trials_analyzed), 'trials', '_PLV.mat'], 'Coh_sec', 'Coh_sec_pct', 'dP_sec', 'dP_sec_pct')

end
