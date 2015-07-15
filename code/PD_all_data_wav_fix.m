function PD_all_data_wav_fix(subjects_mat, freqs, no_cycles, bands)

load(subjects_mat)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, 200);
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    save_name = subj_name;
    
else

    save_name = sprintf('%s_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', subj_name, freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));

end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    no_secs = min(basetime + 1500, length(PD_dec)/sampling_freq);
    
    t = (1:no_secs*sampling_freq)/sampling_freq;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
        load([subj_name, '_wav_laser_artifacts.mat'])
       
        base_index = logical(laser_periods(:, 1));
        
    else
        
        base_index = t <= basetime;
        
    end
    
    clear Spec Spec_norm Spec_pct BP BP_norm BP_pct
    
    Spec = load([save_name, '_wt.mat'], 'Spec');
    
    Spec = abs(Spec.Spec);
    
    [Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(no_secs*sampling_freq, no_freqs, 2));
    
    [BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(no_secs*sampling_freq, no_bands, 2));
    
    for ch = 1:2
        
        %% Baseline normalize.
        
        baseline_mean = mean(abs(Spec(base_index, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Normalize by total power.
        
        Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(sqrt(sum(abs(Spec(:, :, ch)).^2, 2)), 1, no_freqs);
        
        %% Baseline normalize percent of total power.
        
        baseline_mean = mean(abs(Spec_norm(base_index, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Band power.
        
        for b = 1:no_bands
            
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
            
            BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
            
            BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2);
            
            BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
            
        end
        
        BP_pct_baseline_mean = mean(BP_pct(base_index, :, ch));
        
        BP_pct(:, :, ch) = BP_pct(:, :, ch) - ones(size(BP_pct(:, :, ch)))*diag(BP_pct_baseline_mean);
        
        BP_norm_pct_baseline_mean = mean(BP_norm_pct(base_index, :, ch));
        
        BP_norm_pct(:, :, ch) = BP_norm_pct(:, :, ch) - ones(size(BP_norm_pct(:, :, ch)))*diag(BP_norm_pct_baseline_mean);
        
    end
    
    % save([subj_name, '_wt.mat'], '-v7.3', 'Spec', 'Spec_norm', 'Spec_pct', 'Spec_norm_pct', 't', 'basetime', 'freqs', 'sampling_freq')
    
    save([save_name, '_wt.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec', '-v7.3')
    save([save_name, '_wt_norm.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm', '-v7.3')
    save([save_name, '_wt_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_pct', '-v7.3')
    save([save_name, '_wt_norm_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm_pct', '-v7.3')
    
    save([subj_name, '_wt_BP.mat'], '-v7.3', 'BP', 'BP_norm', 'BP_pct', 'BP_norm_pct', 't', 'basetime', 'bands', 'sampling_freq')
    
end