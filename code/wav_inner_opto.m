function wav_inner_opto(folder, prefix, freqs, no_cycles, bands)

subj_name = [folder,'/',prefix];

load([subj_name,'_all_channel_data_dec.mat'])

no_freqs = length(freqs); no_bands = size(bands, 1);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', sampling_freq);

segment_length = sampling_freq;

t = (1:length(PD_dec))/sampling_freq;

clear Spec Spec_norm Spec_pct BP BP_norm BP_pct

[Spec, Spec_norm] = deal(nan(length(PD_dec), no_freqs, 2));

[BP, BP_norm] = deal(nan(length(PD_dec), no_bands, 2));

for ch = 1:2
    
    data = PD_dec(:, ch);
    
    data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))];
    
    for f = 1:no_freqs
        
        conv_prod = conv(data_reflected, wavelets(f,:), 'same');
        
        Spec(:, f, ch) = abs(conv_prod((segment_length + 1):(end - segment_length)));
        
    end
    
    %% Normalize by total power.
    
    Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(sqrt(sum(Spec(:, :, ch).^2, 2)), 1, no_freqs);
    
end

save([subj_name, '_wt.mat'], 'sampling_freq', 't', 'freqs', 'Spec', 'Spec_norm', '-v7.3')
    
%% Band power.

for ch = 1:2

    for b = 1:no_bands
        
        band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
        
        BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
        
    end
    
end

save([subj_name, '_wt_BP.mat'], 'sampling_freq', 't', 'bands', 'BP', 'BP_norm', '-v7.3')
