function wav_inner(folder, prefix, basetime, freqs, no_cycles, bands)

subj_name = [folder,'/',prefix];

load([subj_name,'_all_channel_data_dec.mat'])

no_freqs = length(freqs); no_bands = size(bands, 1);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', sampling_freq);

segment_length = sampling_freq;

t = (1:length(PD_dec))/sampling_freq - basetime;

clear Spec Spec_norm Spec_pct BP BP_norm BP_pct

[Spec, Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(length(PD_dec), no_freqs, 2));

[BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(length(PD_dec), no_bands, 2));

for ch = 1:2
    
    data = PD_dec(:, ch);
    
    data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))];
    
    for f = 1:no_freqs
        
        conv_prod = conv(data_reflected, wavelets(f,:), 'same');
        
        Spec(:, f, ch) = abs(conv_prod((segment_length + 1):(end - segment_length)));
        
    end
    
    %% Baseline normalize.
    
    baseline_mean = mean(abs(Spec(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
    
    Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
    
    %% Normalize by total power.
    
    Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(sqrt(sum(Spec(:, :, ch).^2, 2)), 1, no_freqs);
    
    %% Baseline normalize percent of total power.
    
    baseline_mean = mean(abs(Spec_norm(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
    
    Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
    
end

save([subj_name, '_wt.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec', 'Spec_norm', 'Spec_pct', 'Spec_norm_pct', '-v7.3')
    
%% Band power.

for ch = 1:2

    for b = 1:no_bands
        
        band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
        
        BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
        
        BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
        
        BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
        
    end
    
end

save([subj_name, '_wt_BP.mat'], 'sampling_freq', 't', 'basetime', 'bands', 'BP', 'BP_norm', 'BP_pct', 'BP_norm_pct', '-v7.3')
