function PD_stats_wav_BP(subjects_mat, window_secs, outlier_lims)

freqs = 1:200; no_cycles = linspace(3, 21, 200);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

load(subjects_mat)

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    % load([subj_name, '_', num2str(step_secs), 's_steps_', num2str(window_secs), 's_windows_wav_BP_artifacts.mat'])
    
    load([subj_name, '_wt.mat'], 'sampling_freq')
    
    All_data = load([subj_name, '_wt_BP.mat']);
    
    load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
    
    load([subj_name, '_peaks.mat'])
    
    [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
    
    [~, BP_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
    
    t = All_data.t;
    
    no_windows = floor(t(end)/window_secs);
    
    [t_win, median_pow, mean_pow] = deal(nan(no_windows, 1));
    
    [p_less, p_greater] = deal(nan(no_windows, 2, no_bands, no_norms));
    
    for n = 1:no_norms
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        BP_data(logical(outlier_nans)) = nan;
        
        BP_data(logical(BP_nans)) = nan;
        
        for ch = 1:2
            
            %% Comparing Band Power against Baseline, for Each Window.
            
            for b = 1:no_bands
                
                baseline_data = BP_data(t < 0, b, ch);
                
                baseline_data(isnan(baseline_data)) = [];
                
                for w = 1:no_windows
                    
                    win_indices = t >= (w - 1)*window_secs & t < w*window_secs;
                    
                    t_win(w) = median(t(win_indices));
                    
                    window_data = BP_data(win_indices, b, ch);
                    
                    window_data(isnan(window_data)) = [];
                    
                    if ~isempty(baseline_data) && ~isempty(window_data)
                        
                        median_pow(w) = median(window_data);
                        
                        mean_pow(w) = mean(window_data);
                        
                        p_less(w, ch, b, n) = ranksum(baseline_data, window_data, 'tail', 'right');
                        
                        p_greater(w, ch, b, n) = ranksum(baseline_data, window_data, 'tail', 'left');
                        
                    end
                    
                end
                
            end
            
        end
        
        save([subj_name, '_wt_BP_', num2str(window_secs), 's_stats.mat'], 't_win', 'bands', 'norms', 'median_pow', 'mean_pow', 'p_less', 'p_greater')
        
    end

end

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