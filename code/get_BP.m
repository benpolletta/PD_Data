function BP = get_BP(subj_name, outlier_lim, norm)

load([subj_name, '_wt.mat'], 'sampling_freq', 'freqs')

BP_data = load([subj_name, '_wt_BP.mat']);

bands = BP_data.bands;

BP = getfield(BP_data, ['BP', norm]);

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    [~, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, linspace(3, 21, 200), bands);
    
    laser_nans = repmat(laser_nans, [1 1 2]);
    
    BP(logical(laser_nans)) = nan;
    
end

if ~isempty(outlier_lim) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat']))
    
    load([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat'])
    
    [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, linspace(3, 21, 200), bands);
    
    outlier_nans = repmat(outlier_nans, [1 1 2]);
    
    BP(logical(outlier_nans)) = nan;
    
end

if ~isempty(dir([subj_name, '_peaks.mat']))
    
    load([subj_name, '_peaks.mat'])
    
    [~, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, linspace(3, 21, 200), bands);
    
    BP(logical(spike_nans)) = nan;
    
end

end