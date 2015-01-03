function PD_beta_blocks_rel_infusion(subject_mat, sd_lim)
    
load(subject_mat)

freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

[BP_high_dps, BP_high_cum_dps] = deal(nan(length(folders), 6, no_norms, 2, 2));

for fo = 2:2 %length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];

    load([subj_name, '_wt_BP.mat'], 'sampling_freq')
    
    base_index = basetimes(fo)*sampling_freq;
    
    load([subj_name, '_peaks.mat'])
    
    [~, BP_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, linspace(3, 21, 200), bands);
    
    for n = 1:no_norms
        
        BP_data = load([subj_name,'_wt_BP.mat'], ['BP', norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', norms{n}]);
        
        BP_data(logical(BP_nans)) = nan;
        
        t = (1:size(BP_data, 1))/sampling_freq;
        
        BP_high = nan(size(BP_data));
        
        pd_limits = [1 min(length(t), base_index); min(length(t), base_index + 1) min(length(t), base_index + 1500*sampling_freq)];
        
        pd_lengths = diff(pd_limits, [], 2) + 1;
        
        for ch = 1:2
            
            for b = 1:size(BP_high, 2)
                
                high_cutoff = nanmean(BP_data(1:pd_limits(1, 2), b, ch)) + sd_lim*nanstd(BP_data(1:pd_limits(1, 2), b, ch));
                
                BP_high(:, b, ch) = BP_data(:, b, ch) >= high_cutoff;
                
            end
            
        end
        
        BP_high_cum = BP_high;
        
        BP_high_cum(cumsum(BP_high, 2) > 1) = 0;
        
        for pd = 1:size(pd_limits,1)
            
            for ch = 1:2
                
                BP_high_dps(fo, :, ch, pd, n) = sum(BP_high(pd_limits(pd,1):pd_limits(pd,2), :, ch))/pd_lengths(pd);
                    
                BP_high_cum_dps(fo, :, ch, pd, n) = sum(BP_high_cum(pd_limits(pd,1):pd_limits(pd,2), :, ch))/pd_lengths(pd);
                
            end
            
        end
        
        save([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
        
    end
          
end

save([subject_mat(1:(end - length('_subjects.mat'))), '_', num2str(sd_lim), 'sd_BP_high_dps.mat'], 'BP_high_dps', 'BP_high_cum_dps')

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