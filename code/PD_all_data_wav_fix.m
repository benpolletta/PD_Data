function PD_all_data_wav_fix(subjects_mat)

load(subjects_mat)

% sampling_freq = 1000;

freqs = 1:200; % no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

% no_cycles = linspace(3, 21, no_freqs);
% 
% wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', sampling_freq);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    base_index = basetime*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    no_secs = min(basetime + 1500, length(PD_dec)/sampling_freq);
    
    t = (1:no_secs*sampling_freq)/sampling_freq;
    
    clear Spec Spec_norm Spec_pct BP BP_norm BP_pct
    
    Spec = load([subj_name, '_wt.mat'], 'Spec');
    
    Spec = abs(Spec.Spec);
    
    [Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(no_secs*sampling_freq, no_freqs, 2));
    
    [BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(no_secs*sampling_freq, no_bands, 2));
    
    for ch = 1:2
        
        %% Baseline normalize.
        
        baseline_mean = mean(abs(Spec(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Normalize by total power.
        
        Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(sqrt(sum(abs(Spec(:, :, ch)).^2, 2)), 1, no_freqs);
        
        %% Baseline normalize percent of total power.
        
        baseline_mean = mean(abs(Spec_norm(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Band power.
        
        for b = 1:no_bands
            
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
            
            BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
            
            BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
            
            BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
            
        end
        
    end
    
    save([subj_name, '_wt.mat'], '-v7.3', 'Spec', 'Spec_norm', 'Spec_pct', 'Spec_norm_pct', 't', 'basetime', 'freqs', 'sampling_freq')
    
    save([subj_name, '_wt_BP.mat'], '-v7.3', 'BP', 'BP_norm', 'BP_pct', 'BP_norm_pct', 't', 'basetime', 'bands', 'sampling_freq')
    
end