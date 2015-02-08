function [BP, Spec] = get_BP(subj_name, outlier_lim, norm, freqs, no_cycles, bands)

load([subj_name, '_wt.mat'], 'sampling_freq') % , 'freqs')

BP_data = load([subj_name, '_wt_BP.mat'], ['BP', norm]);

BP = getfield(BP_data, ['BP', norm]);

Spec_data = load([subj_name, '_wt', norm, '.mat'], ['Spec', norm]);

Spec = getfield(Spec_data, ['Spec', norm]);

% bands = BP_data.bands;

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    [laser_nans_wt, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, no_cycles, bands);
    
    laser_nans = repmat(laser_nans, [1 1 2]);
    
    BP(logical(laser_nans)) = nan;
    
    laser_nans_wt = repmat(laser_nans_wt, [1 1 2]);
    
    Spec(logical(laser_nans)) = nan;
    
end

if ~isempty(outlier_lim) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat']))
    
    load([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat'])
    
    [outlier_nans_wt, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
    
    outlier_nans = repmat(outlier_nans, [1 1 2]);
    
    BP(logical(outlier_nans)) = nan;
    
    outlier_nans_wt = repmat(outlier_nans_wt, [1 1 2]);
    
    Spec(logical(outlier_nans_wt)) = nan;
    
end

if ~isempty(dir([subj_name, '_peaks.mat']))
    
    load([subj_name, '_peaks.mat'])
    
    [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
    
    BP(logical(spike_nans)) = nan;
    
    Spec(logical(spike_nans_wt)) = nan;
    
end

end