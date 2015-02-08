function PD_beta_blocks_rel_infusion(subject_mat, sd_lim, outlier_lims, freqs, no_cycles, bands)

% Computes zero-one vectors indicating time points at which band power is
% above sd_lim standard deviations from the (pre-infusion) mean.
%
% INPUTS:
% subjects_mat - name of .mat file containing list of folders.
% sd_lim - standard deviation cutoff above which the power at a given
%  timepoint is deemed "high"; typical value is 2.
% outlier_lims - vector of values, same length as 'folders' cell in
%   subject_mat, which indicates the outlier cutoff used to remove
%   artifacts from each recording; if this is left empty, a default vector
%   included in subject_mat is used.
% freqs - a vector of center frequencies for the wavelet spectrogram;
%  default is 1:200.
% no_cycles - a vector of numbers of cycles each wavelet contains
%  (effectively gives the bandwidth of the spectrogram at each frequency);
%  default is linspace(3, 7, 200).
% bands - a matrix of band limits, over which spectrogram power will be
%  summed; matrix is n by 2, with n the number of bands, the first
%  column containing the start frequency of each band, and the second column
%  containing the end frequency of each band; default is [1 4;4 8;8 30;30
%  100;120 180;0 200].

if isempty(freqs) && isempty(no_cycles) && isempty(bands) % Defaults.
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = '';
    
else % Suffix used to save non-default output.
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
end
    
load(subject_mat)

no_bands = size(bands, 1);

no_norms = length(norms);

no_pds = length(pd_labels);

no_chans = length(chan_labels);

% [BP_high_dps, BP_high_cum_dps] = deal(nan(length(folders), no_bands, no_chans, no_pds, no_norms));

for fo = 1:length(folders)
    
    clear pd_indices base_index
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    no_bands = size(bands, 1);

    load([subj_name, BP_suffix, '_wt_BP.mat'])
    
    base_index = basetimes(fo)*sampling_freq;
        
    t = ((1:size(BP, 1)) - base_index)/sampling_freq;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat'])) % For optogenetic data.
    
        load([subj_name, '_wav_laser_artifacts.mat'])
        
        [~, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, no_cycles, bands); 
        % Transform indicators of laser on/off periods to matrices
        % containing nans where spectrograms/band power should be blocked
        % out.
        
        laser_nans = repmat(laser_nans, [1 1 2]); % Same matrix for each channel.
        
        pd_indices = laser_periods; % Matrix containing periods for analysis.
        
        base_index = logical(ones(size(BP, 1), 1)); % pd_indices(:, 1);
        % Index of timepoints to use when calculating high power cutoff;
        % for optogenetics data, use entire data span.
    
    else
        
        laser_nans = [];
        
        if no_pds == 2 % For carbachol infusion.
        
            pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0; % Matrix containing periods for analysis.
            
            base_index = pd_indices(:, 1);
            % Index of timepoints to use when calculating high power cutoff;
            % for carbachol data, use pre-infusion data.
            
        elseif no_pds == 1
           
            pd_indices = ones(size(t)); % Matrix containing periods for analysis.
            
            base_index = logical(ones(size(BP, 1), 1));
            % Index of timepoints to use when calculating high power cutoff;
            % for 6OHDA data, use entire data span.
            
        end
        
    end
    
    pd_indices = logical(pd_indices); base_index = logical(base_index);
    
    if ~isempty(outlier_lims) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])) % Removing outliers.
    
        load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
        
        [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
        
        outlier_nans = repmat(outlier_nans, [1 1 2]);
        
    else
        
        outlier_nans = [];
        
    end
    
    if ~isempty(dir([subj_name, '_peaks.mat'])) % Removing peaks.
        
        load([subj_name, '_peaks.mat'])
        
        [~, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
        
    else
        
        spike_nans = [];
        
    end
    
    for n = 1:no_norms
        
        BP_data = load([subj_name, BP_suffix, '_wt_BP.mat'], ['BP', norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', norms{n}]);
        
        if ~isempty(laser_nans)
            
            BP_data(logical(laser_nans)) = nan;
            
        end
        
        if ~isempty(outlier_nans)
            
            BP_data(logical(outlier_nans)) = nan;
            
        end
        
        if ~isempty(spike_nans)
            
            BP_data(logical(spike_nans)) = nan;
            
        end
        
        BP_high = nan(size(BP_data));
        
        % base_limits = [1 min(length(t), floor(base_index))]; % ; min(length(t) - 1, floor(base_index) + 1) length(t)];
        
        for ch = 1:no_chans
            
            for b = 1:size(BP_high, 2)
                
                high_cutoff = nanmean(BP_data(base_index, b, ch)) + sd_lim*nanstd(BP_data(base_index, b, ch));
                
                BP_high(:, b, ch) = BP_data(:, b, ch) >= high_cutoff;
                
            end
            
        end
        
        BP_high_cum = BP_high;
        
        BP_high_cum(cumsum(BP_high, 2) > 1) = 0;
        
        save([subj_name, BP_suffix, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
        
        % for pd = 1:no_pds
        % 
        %     for ch = 1:no_chans
        % 
        %         BP_high_dps(fo, :, ch, pd, n) = nansum(BP_high(pd_indices(:, pd), :, ch))/sum(pd_indices(:, pd));
        % 
        %         BP_high_cum_dps(fo, :, ch, pd, n) = nansum(BP_high_cum(pd_indices(:, pd), :, ch))/sum(pd_indices(:, pd));
        % 
        %     end
        % 
        % end
        
    end
          
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_', num2str(sd_lim), 'sd_BP_high_dps.mat'], 'BP_high_dps', 'BP_high_cum_dps')

end

function [wav_nans, BP_nans] = indicator_to_nans(indicator, sampling_freq, freqs, no_cycles, bands)

[no_dps, no_channels] = size(indicator);

wav_nans = nan(no_dps, length(freqs), no_channels);

BP_nans = nan(no_dps, size(bands, 1), no_channels);

for ch = 1:no_channels
    
    wav_nans_temp = abs(wavelet_spectrogram(indicator(:, ch), sampling_freq, freqs, no_cycles, 0, ''));
    
    wav_nans_temp = wav_nans_temp*diag(1./max(wav_nans_temp));
    
    wav_nans(:, :, ch) = wav_nans_temp > .01;
    
    for b = 1:size(bands, 1)
       
        band_freqs = freqs >= bands(b, 1) & freqs <= bands(b, 2);
        
        BP_nans(:, b, ch) = sum(wav_nans(:, band_freqs, ch), 2);
        
    end
    
    BP_nans(BP_nans > 0) = 1;

end
    
end