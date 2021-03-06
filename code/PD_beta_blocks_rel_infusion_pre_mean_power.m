function PD_beta_blocks_rel_infusion_pre_mean_power(subject_mat, epoch_secs, pd_handle, norm, freqs, no_cycles, bands)

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

no_folders = length(folders);

load([folders{1}, '/', prefixes{1}, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

BP_pre_mean = nan(no_folders, no_bands, 2);

load([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

pre_mean_name = [subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', norm, '_power_pre'];

pre_mean_fid = nan(no_chans, 1);

for ch = 1:no_chans

    pre_mean_fid(ch) = fopen([pre_mean_name, '_ch', num2str(ch), '.txt'], 'w');

    fprintf(pre_mean_fid(ch), make_format(no_bands + 1, 's'), 'Recording', short_band_labels{:});
    
end

pre_mean_format = make_format(no_bands, 'f');

pre_mean_format = ['%s\t', pre_mean_format];

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    if strcmp(norm, '_peaks')
        
        load([subj_name, '_peaks.mat'])
        
        BP = repmat(permute(Spike_indicator, [1 2 3]), [1 no_bands 1]);
        
    else
    
        load([subj_name, BP_suffix, '_wt_BP.mat'])
        
	eval(['BP = BP', norm, ';'])
        
    end
    
    t = (1:size(BP, 1))/sampling_freq - basetimes(fo); % t = t/60;
    
    clear pd_indices
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), 1);
        
    end
    
    pd_indices = logical(pd_indices);
    
    if ~strcmp(norm, '_spikes')
        
        if ~isempty(outlier_lims) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat']))
            
            load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
            
            [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
            
            outlier_nans = repmat(outlier_nans, [1 1 2]);
            
            BP(logical(outlier_nans)) = nan;
            
        end
        
        if ~isempty(dir([subj_name, '_peaks.mat']))
            
            load([subj_name, '_peaks.mat'])
            
            [~, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
            
            BP(logical(spike_nans)) = nan;
            
        end
        
    end
        
    BP_pre_mean(fo, :, :) = nanmean(BP(pd_indices(:, 1), :, :));
    
    for ch = 1:no_chans
        
        fprintf(pre_mean_fid(ch), pre_mean_format, folder, BP_pre_mean(fo, :, ch));
        
    end
       
end

for ch = 1:no_chans, fclose(pre_mean_fid(ch)); end

save([pre_mean_name, '.mat'], 'BP_pre_mean')

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
