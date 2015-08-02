function [wav_nans, BP_nans] = indicator_to_nans(indicator, sampling_freq, freqs, no_cycles, bands)

[no_dps, no_channels] = size(indicator);

wav_nans = nan(no_dps, length(freqs), no_channels);

BP_nans = nan(no_dps, size(bands, 1), no_channels);

for ch = 1:no_channels
    
    wav_nans_temp = abs_wavelet_spectrogram(indicator(:, ch), sampling_freq, freqs, no_cycles, 0, '');
    
    wav_nans_temp = wav_nans_temp*diag(1./max(wav_nans_temp));
    
    wav_nans(:, :, ch) = wav_nans_temp > .01;
    
    for b = 1:size(bands, 1)
       
        band_freqs = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP_nans(:, b, ch) = sum(wav_nans(:, band_freqs, ch), 2);
        
    end
    
    BP_nans(BP_nans > 0) = 1;

end
    
end

% FUNCTION FROM INSIDE PD_REL_INFUSION_PLOT_SPECTROGRAM, FOR LATER CHECKING
% WITH ABOVE (8/1/2015).
% 
% function [wav_nans, BP_nans] = indicator_to_nans(indicator, sampling_freq, freqs, no_cycles, bands)
% 
% [no_dps, no_channels] = size(indicator);
% 
% wav_nans = nan(no_dps, length(freqs), no_channels);
% 
% BP_nans = nan(no_dps, size(bands, 1), no_channels);
% 
% for ch = 1:no_channels
%     
%     wav_nans_temp = abs(wavelet_spectrogram(indicator(:, ch), sampling_freq, freqs, no_cycles, 0, ''));
%     
%     wav_nans_temp = wav_nans_temp*diag(1./max(wav_nans_temp));
%     
%     wav_nans(:, :, ch) = wav_nans_temp > .01;
%     
%     for b = 1:size(bands, 1)
%        
%         band_freqs = freqs >= bands(b, 1) & freqs <= bands(b, 2);
%         
%         BP_nans(:, b, ch) = sum(wav_nans(:, band_freqs, ch), 2);
%         
%     end
%     
%     BP_nans(BP_nans > 0) = 1;
% 
% end
%     
% end