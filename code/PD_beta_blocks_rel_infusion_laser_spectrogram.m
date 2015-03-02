function PD_beta_blocks_rel_infusion_laser_spectrogram(subject_mat, norm, no_trials_analyzed, band_index, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
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
    
secs_per_trial = 5;

timepoints_per_sec = 20;

timepoints_per_trial = secs_per_trial*timepoints_per_sec;

dps_per_timepoint = sampling_freq/timepoints_per_sec; % no_timepoints = max_no_trials*secs_per_trial*timepoints_per_sec;

WT_spec = nan(sum(band_indices{band_index}), timepoints_per_trial, no_pds, no_chans, no_folders, max_no_trials);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    [~, Spec_data] = get_BP(subj_name, outlier_lims(fo), norm, freqs, no_cycles, bands);
    
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
                
                for tp = 1:(timepoints_per_sec*secs_per_trial)
                    
                    tp_start = trial_start + (tp - 1)*dps_per_timepoint + 1;
                    
                    tp_end = trial_start + tp*dps_per_timepoint;
                    
                    WT_spec(:, tp, pd, ch, fo, tr) =...
                        nanmean(Spec_data(tp_start:tp_end, band_indices{band_index}, ch))';
                    
                end
                
            end
            
        end
        
    end
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_', short_band_labels{band_index},...
    '_', num2str(no_trials_analyzed), 'trials', norm, '_', num2str(timepoints_per_sec), 'Hz_spectrogram.mat'], 'WT_spec')

plot_spectrogram(WT_spec, no_chans)

save_as_pdf(gcf, [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_', short_band_labels{band_index},...
    '_', num2str(no_trials_analyzed), 'trials', norm, '_', num2str(timepoints_per_sec), 'Hz_spectrogram'])

end


function plot_spectrogram(WT_spec, no_chans)

mean_WT_spec = nanmean(WT_spec, 6);

% for fo = 1:size(WT_spec, 4)
%    
%     folder_spec = mean_WT_spec(
%     
% end

mean_mean_WT_spec = nanmean(mean_WT_spec, 5);

figure

for ch = 1:no_chans
    
    subplot(no_chans, 1, ch)
    
    chan_spec = mean_mean_WT_spec(:, :, :, ch);
    
    chan_spec = reshape(chan_spec, size(chan_spec, 1), size(chan_spec, 2)*size(chan_spec, 3));
    
    imagesc(chan_spec)
    
    axis xy
    
end

end
