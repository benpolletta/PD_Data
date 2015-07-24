function PD_pct_fix(subjects_mat, freqs, no_cycles, bands)

% Fix percent normalization, especially for laser data (7/15/15).

load(subjects_mat)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3, 21, 200); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = '';
    
else

    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
    BP_suffix = sprintf('%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));

end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    % no_secs = min(basetime + 1500, length(PD_dec)/sampling_freq);
    
    t = (1:no_secs*sampling_freq)/sampling_freq;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
        load([subj_name, '_wav_laser_artifacts.mat'])
       
        base_index = laser_periods(:, 1);

        base_trials = index_to_blocks(base_index);
        
        base_end = base_trials(10, 2);
        
        base_index((base_end + 1):end) = 0;
        
        base_index = logical(base_index);
        
    else
        
        base_index = t <= basetime;
        
    end
    
    clear Spec Spec_pct BP BP_pct
    
    [~, Spec] = get_BP(subj_name, outlier_lims(fo), '', in_freqs, in_no_cycles, in_bands);

    % Spec = load([subj_name, BP_suffix, '_wt.mat'], 'Spec');
    % 
    % Spec = abs(Spec.Spec);
    
    Spec_pct = nan(size(Spec, 1), length(freqs), 2);
    
    BP_pct = nan(size(Spec, 1), size(bands, 1), 2);
    
    for ch = 1:2
        
        %% Baseline normalize.
        
        baseline_mean = nanmean(abs(Spec(base_index, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
        
        Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
        
        %% Band power.
        
        for b = 1:size(bands, 1)
            
            band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
            
            BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2);
            
        end
        
        BP_pct_baseline_mean = mean(BP_pct(base_index, :, ch));
        
        BP_pct(:, :, ch) = BP_pct(:, :, ch) - ones(size(BP_pct(:, :, ch)))*diag(BP_pct_baseline_mean);
        
    end
    
    save([subj_name, BP_suffix, '_wt_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'no_cycles', 'bands', 'Spec_pct', '-v7.3')
    
    save([subj_name, BP_suffix, '_wt_BP.mat'], '-v7.3', '-append', 'BP_pct')
    
end