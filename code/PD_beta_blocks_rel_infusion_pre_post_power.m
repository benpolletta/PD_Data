function PD_beta_blocks_rel_infusion_pre_post_power(subject_mat, epoch_secs, pd_handle)
    
close('all')

load(subject_mat)

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

BP_sec = nan(epoch_secs, no_folders, 2, no_bands, 2);

load([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_wt_BP.mat'])
    
    t = (1:size(BP, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    clear pd_indices
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
        load([subj_name, '_wav_laser_artifacts.mat'])
        
        [~, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, linspace(3, 21, 200), bands);
        
        laser_nans = repmat(laser_nans, [1 1 2]);
        
        BP(logical(laser_nans)) = nan;
        
        pd_indices = laser_periods;
    
    else
        
        if no_pds == 2
        
            pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
            
        elseif no_pds == 1
           
            pd_indices = ones(length(t), 1);
            
        end
        
    end
    
    pd_indices = logical(pd_indices);
    
    if ~isempty(outlier_lims) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat']))
    
        load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
        
        [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, linspace(3, 21, 200), bands);
        
        outlier_nans = repmat(outlier_nans, [1 1 2]);
        
        BP(logical(outlier_nans)) = nan;
        
    end
    
    if ~isempty(dir([subj_name, '_peaks.mat']))
        
        load([subj_name, '_peaks.mat'])
        
        [~, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, linspace(3, 21, 200), bands);
        
        BP(logical(spike_nans)) = nan;
        
    end
    
    % beta_blocks_find = nan(size(BP_high_cum));
    
    for b = 1:no_bands
        
        figure(b)
        
        for ch = 1:no_chans
            
            handle = nan(no_pds, 1);
            
            for pd = 1:no_pds
                
                bp_max_start = All_bp_max_start(fo, ch, b, pd);
                
                bp_max_end = All_bp_max_end(fo, ch, b, pd);
                
                BP_plot = nanconv(BP(:, b, ch), ones(60*sampling_freq, 1)/(60*sampling_freq), 'nanout');
                
                subplot(no_folders, 2, (fo - 1)*2 + ch)
                
                plot(t/60, BP_plot)
                
                axis tight
                
                hold on
                
                handle(pd) = plot(t(bp_max_start:bp_max_end)/60, BP_plot(bp_max_start:bp_max_end), pd_colors{pd}, 'LineWidth', 2);
                
                for sec = 1:epoch_secs
                   
                    sec_start = max(bp_max_start + (sec - 1)*sampling_freq + 1, 1);
                    
                    sec_end = min(max(bp_max_start + sec*sampling_freq, 1), find(pd_indices(:, pd) == 1, 1, 'last'));
                    
                    BP_sec(sec, fo, pd, b, ch) = nanmean(BP(sec_start:sec_end, b, ch));
                    
                end
                
            end
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
                
            elseif fo == no_folders
                
                xlabel('Time Rel. Infusion (Min.)')
                
            end
            
            if ch == 1
                
                ylabel({folder; 'High Power Density per Min.'})
                
                if fo == 1
                    
                    for pd = 1:no_pds
                        
                        ledge{pd} = ['Peak ', num2str(epoch_secs/60), ' Min., ', pd_labels{pd}];
                        
                    end
                    
                    legend(handle, ledge)
                    
                end
                
            end
            
            plot([0; 0], [min(BP_plot); max(BP_plot)], 'k', 'LineWidth', 1)
            
        end
        
    end
    
end

for b = 1:no_bands
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_min_', short_band_labels{b}, pd_handle, '_power'])
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '_power.mat'], 'BP_sec')

end

function [wav_nans, BP_nans] = indicator_to_nans(indicator, sampling_freq, freqs, no_cycles, bands)

[no_dps, no_channels] = size(indicator);

wav_nans = nan(no_dps, length(freqs), no_channels);

BP_nans = nan(no_dps, size(bands, 1), no_channels);

for ch = 1:no_channels
    
    wav_nans_temp = abs(wavelet_spectrogram(indicator(:, ch), sampling_freq, freqs, no_cycles, 0, ''));
    
    wav_nans_temp = wav_nans_temp*diag(1./max(wav_nans_temp));
    
    wav_nans(:, :, ch) = wav_nans_temp > .01;
    
    for b = 1:size(bands, 1)
       
        band_freqs = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP_nans(:, b, ch) = sum(wav_nans(:, band_freqs, ch), 2);
        
    end
    
    BP_nans(BP_nans > 0) = 1;

end
    
end