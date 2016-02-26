function [BP, Spec] = get_BP(subj_name, peak_suffix, outlier_lim, norm, freqs, no_cycles, bands)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 7, 200);
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end

load([subj_name, BP_suffix, '_wt.mat'], 'sampling_freq') % , 'freqs')

BP_data = load([subj_name, BP_suffix, '_wt_BP.mat'], ['BP', norm]);

BP = getfield(BP_data, ['BP', norm]);

Spec_data = load([subj_name, BP_suffix, '_wt', norm, '.mat'], ['Spec', norm]);

Spec = getfield(Spec_data, ['Spec', norm]);

% bands = BP_data.bands;

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    [laser_nans_wt, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, no_cycles, bands);
    
    laser_nans = repmat(laser_nans, [1 1 2]);
    
    BP(logical(laser_nans)) = nan;
    
    laser_nans_wt = repmat(laser_nans_wt, [1 1 2]);
    
    Spec(logical(laser_nans_wt)) = nan;
    
end

if ~isempty(outlier_lim) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat']))
    
    load([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat'])
    
    [outlier_nans_wt, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
    
    outlier_nans = repmat(outlier_nans, [1 1 2]);
    
    BP(logical(outlier_nans)) = nan;
    
    outlier_nans_wt = repmat(outlier_nans_wt, [1 1 2]);
    
    Spec(logical(outlier_nans_wt)) = nan;
    
end

if isempty(peak_suffix)
    
    if ~isempty(dir([subj_name, '_peaks.mat']))
        
        load([subj_name, '_peaks.mat'])
        
        [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
        
        BP(logical(spike_nans)) = nan;
        
        Spec(logical(spike_nans_wt)) = nan;
        
    elseif ~isempty(dir([subj_name, '_chan1_artifacts.mat'])) || ~isempty(dir([subj_name, '_chan2_artifacts.mat']))
        
        for ch = 1:2
            
            load([subj_name, '_chan', num2str(ch), '_artifacts.mat'])
            
            Spike_indicator(:, ch) = peak_indicator;
            
        end
        
        [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
        
        BP(logical(spike_nans)) = nan;
        
        Spec(logical(spike_nans_wt)) = nan;
        
    end

elseif strcmp(peak_suffix, '_kmeans')
        
    if ~isempty(dir([subj_name, '_chan1_artifacts.mat'])) || ~isempty(dir([subj_name, '_chan2_artifacts.mat']))
        
        for ch = 1:2
            
            load([subj_name, '_chan', num2str(ch), '_artifacts.mat'])
            
            Spike_indicator(:, ch) = peak_indicator;
            
        end
        
        [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
        
        BP(logical(spike_nans)) = nan;
        
        Spec(logical(spike_nans_wt)) = nan;
        
    elseif ~isempty(dir([subj_name, '_peaks.mat']))
        
        load([subj_name, '_peaks.mat'])
        
        [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
        
        BP(logical(spike_nans)) = nan;
        
        Spec(logical(spike_nans_wt)) = nan;
        
    end
    
end