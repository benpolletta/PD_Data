function PD_pmtm_normalize(subjects_mat, epoch_secs)

load(subjects_mat)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 500]; no_bands = size(bands, 1);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_', num2str(epoch_secs),'s_epoch_pmtm.mat'], 't', 'freqs', 'Spec', 'Spec_norm')
        
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.mat'])
    
    Spec(artifact_indicator, :, :) = nan;
    
    Spec_norm(artifact_indicator, :, :) = nan;
    
    [no_epochs, no_freqs, ~] = size(Spec);
    
    [Spec_pct, Spec_norm_pct] = deal(nan(no_epochs, no_freqs, 2));
        
    [BP_pct, BP_norm_pct] = deal(nan(no_epochs, no_bands, 2));
    
    for ch = 1:2
        
        %% Baseline normalize.
        
        baseline_mean = nanmean(abs(Spec(t < 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Baseline normalize percent of total power.
        
        baseline_mean = nanmean(abs(Spec_norm(t < 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Band power.
        
        for b = 1:no_bands
                
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
            
            BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2);
            
        end
        
    end
        
    save([subj_name, '_', num2str(epoch_secs),'s_epoch_pmtm.mat'], '-v7.3', '-append', 'Spec_pct', 'Spec_norm_pct', 'BP_pct', 'BP_norm_pct')
       
end