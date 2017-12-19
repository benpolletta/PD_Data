function wav_inner(folder, prefix, basetime, freqs, no_cycles, bands, peak_suffix)

subj_name = [folder,'/',prefix];

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, 200);
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    save_name = [subj_name, peak_suffix];
    
else

    save_name = sprintf('%s_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands%s', subj_name,...
        freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1), peak_suffix);

end
    
load([subj_name,'_all_channel_data_dec.mat'])

no_freqs = length(freqs); no_bands = size(bands, 1);

cycle_lengths = no_cycles.*(sampling_freq./freqs);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', max(sampling_freq, max(cycle_lengths)));

segment_length = sampling_freq;

t = (1:length(PD_dec))/sampling_freq - basetime;

clear Spec Spec_norm Spec_pct Spec_norm_pct BP BP_norm BP_pct BP_norm_pct

%% Wavelet spectrogram.

[Spec, Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(length(PD_dec), no_freqs, 2));

for ch = 1:2
    
    data = PD_dec(:, ch);
    
    data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))];
    
    for f = 1:no_freqs
        
        conv_prod = conv(data_reflected, wavelets(f,:), 'same');
        
        Spec(:, f, ch) = conv_prod((segment_length + 1):(end - segment_length));
        
    end
    
    %% Baseline normalize.
    
    baseline_mean = mean(abs(Spec(t <= 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
    
    Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
    
    %% Normalize by total power.
    
    Spec_norm(:, :, ch) = abs(Spec(:, :, ch))./repmat(sqrt(sum(abs(Spec(:, :, ch)).^2, 2)), 1, no_freqs);
    
    %% Baseline normalize percent of total power.
    
    baseline_mean = mean(abs(Spec_norm(t <= 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
    
    Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
    
end

% save([subj_name, '_wt.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec', 'Spec_norm', 'Spec_pct', 'Spec_norm_pct', '-v7.3')

save([save_name, '_wt.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec', '-v7.3')
save([save_name, '_wt_norm.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm', '-v7.3')
save([save_name, '_wt_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_pct', '-v7.3')
save([save_name, '_wt_norm_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm_pct', '-v7.3')
    
%% Band power.

[BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(length(PD_dec), no_bands, 2));

for ch = 1:2

    for b = 1:no_bands
        
        band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
        
        % BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
        
        BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
        
        % BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
        
    end
    
end

save([save_name, '_wt_BP.mat'], 'sampling_freq', 't', 'basetime', 'bands', 'BP', 'BP_norm', '-v7.3') % 'BP_norm', 'BP_pct', 'BP_norm_pct', '-v7.3')
