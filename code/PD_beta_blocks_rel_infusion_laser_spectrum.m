function PD_beta_blocks_rel_infusion_laser_spectrum(subject_mat, peak_suffix, norm, no_trials_analyzed, band_index, freqs, no_cycles, bands)

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
        
        trials{fo, pd}(trials{fo, pd}(1,1) < trials{fo, 1}(1,1), :) = [];
        
        no_trials(fo, pd) = size(trials{fo, pd}, 1);
        
    end
    
end
    
no_trials = min(no_trials, [], 2);

if isempty(no_trials_analyzed)

    max_no_trials = all_dimensions(@max, no_trials);

else
    
    max_no_trials = no_trials_analyzed;
    
end
    
no_secs = max_no_trials*5;

[WT_sec, WT_trial_normed_sec] = deal(nan(sum(band_indices{band_index}), no_secs, no_folders, no_pds, no_chans));
        
[WT_trial_normed_mean, WT_trial_normed_se] = deal(nan(sum(band_indices{band_index}), no_folders, no_pds, no_chans));

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    [~, Spec_data] = get_BP(subj_name, peak_suffix, no_trials, outlier_lims(fo), norm, in_freqs, in_no_cycles, in_bands);
    
    if strcmp(norm, '')
        
        Spec_data = abs(Spec_data);
        
    end
    
    for ch = 1:no_chans
        
        [WT_trials, WT_trials_normed] = deal(nan(min(max_no_trials, no_trials(fo))*5*sampling_freq, length(band_indices{band_index}), no_pds));
        
        for tr = 1:min(max_no_trials, no_trials(fo))
        
            for pd = 1:no_pds
                
                trial_start = trials{fo, pd}(tr, 1);
                
                % for sec = 1:5
                % 
                %     sec_start = trial_start + (sec - 1)*sampling_freq;
                % 
                %     sec_end = trial_start + sec*sampling_freq - 1;
                % 
                %     WT_sec(:, (tr - 1)*5 + sec, fo, pd, ch) = nanmean(Spec_data(sec_start:sec_end, band_indices{band_index}, ch))';
                % 
                % end
                
                trial_end = trial_start + 5*sampling_freq - 1;
                
                WT_trial = Spec_data(trial_start:trial_end, band_indices{band_index}, ch);
                
                WT_trials((tr - 1)*5*sampling_freq + (1:5*sampling_freq), :, pd) = WT_trial;
                
                if pd == 1
                    
                    trial_baseline = nanmean(WT_trial);
                    
                end
                
                WT_trials_normed((tr - 1)*5*sampling_freq + (1:5*sampling_freq), :, pd) = 100*WT_trial./(ones(size(WT_trial))*diag(trial_baseline)) - 100;
                
            end
            
        end
                
        for pd = 1:no_pds
                
            WT_trials(:, :, pd) = nans_to_end(WT_trials(:, :, pd));
            
            WT_trials_normed(:, :, pd) = nans_to_end(WT_trials_normed(:, :, pd));
            
            for sec = 1:min(max_no_trials, no_trials(fo))*5
               
                sec_start = (sec - 1)*sampling_freq + 1;
                
                sec_end = sec*sampling_freq;
                
                WT_sec(:, sec, fo, pd, ch) = nanmean(WT_trials(sec_start:sec_end, :, pd))';
                
                WT_trial_normed_sec(:, sec, fo, pd, ch) = nanmean(WT_trials_normed(sec_start:sec_end, :, pd))';
                
            end
            
            WT_trial_normed_mean(:, fo, pd, ch) = nanmean(WT_trials_normed(:, :, pd))';
            
            WT_trial_normed_se(:, fo, pd, ch) = nanstd(WT_trials_normed(:, :, pd))';
            
        end
        
    end
    
end

if strcmp(norm, '')
    
    save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_', short_band_labels{band_index},...
        '_', num2str(no_trials_analyzed), 'trials', norm, '_spectrum.mat'], 'WT_sec', 'WT_trial_normed_sec', 'WT_trial_normed_mean', 'WT_trial_normed_se')
    
else
    
    save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_', short_band_labels{band_index},...
        '_', num2str(no_trials_analyzed), 'trials', norm, '_spectrum.mat'], 'WT_sec')
    
end

end
