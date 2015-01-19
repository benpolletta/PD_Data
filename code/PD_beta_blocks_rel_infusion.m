function PD_beta_blocks_rel_infusion(subject_mat, sd_lim, outlier_lims)
    
load(subject_mat)

freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

no_norms = length(norms);

no_pds = length(pd_labels);

no_chans = length(chan_labels);

[BP_high_dps, BP_high_cum_dps] = deal(nan(length(folders), no_bands, no_chans, no_pds, no_norms));

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];

    load([subj_name, '_wt_BP.mat'])
    
    base_index = basetimes(fo)*sampling_freq;
        
    t = ((1:size(BP, 1)) - base_index)/sampling_freq;
    
    if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
        load([subj_name, '_wav_laser_artifacts.mat'])
        
        [~, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, linspace(3, 21, 200), bands);
        
        pd_indices = laser_periods;
        
        base_index = logical(ones(size(BP, 1), 1)); % pd_indices(:, 1); %
    
    else
        
        laser_nans = [];
        
        if no_pds == 2
        
            pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
            
            base_index = pd_indices(:, 1);
            
        elseif no_pds == 1
           
            pd_indices = ones(size(t));
            
            base_index = logical(ones(size(BP, 1), 1));
            
        end
        
    end
    
    pd_indices = logical(pd_indices); base_index = logical(base_index);
    
    % if no_pds == 2
    % 
    %     base_index = pd_indices(:, 1);
    % 
    % elseif no_pds
    % 
    %     base_index = logical(ones(size(BP, 1), 1));
    % 
    % end
    
    if ~isempty(outlier_lims) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat']))
    
        load([subj_name, '_wav_BP_', num2str(outlier_lims(fo)), 'sd_outliers.mat'])
        
        [~, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, linspace(3, 21, 200), bands);
        
    else
        
        outlier_nans = [];
        
    end
    
    if ~isempty(dir([subj_name, '_peaks.mat']))
        
        load([subj_name, '_peaks.mat'])
        
        [~, BP_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, linspace(3, 21, 200), bands);
        
    else
        
        BP_nans = [];
        
    end
    
    for n = 1:no_norms
        
        BP_data = load([subj_name,'_wt_BP.mat'], ['BP', norms{n}]);
        
        BP_data = getfield(BP_data, ['BP', norms{n}]);
        
        if ~isempty(laser_nans)
            
            BP_data(logical(laser_nans)) = nan;
            
        end
        
        if ~isempty(outlier_nans)
            
            BP_data(logical(outlier_nans)) = nan;
            
        end
        
        if ~isempty(BP_nans)
            
            BP_data(logical(BP_nans)) = nan;
            
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
        
        save([subj_name, '_', num2str(sd_lim), 'sd_BP', norms{n}, '_high.mat'], 'BP_high', 'BP_high_cum')
        
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

save([subject_mat(1:(end - length('_subjects.mat'))), '_', num2str(sd_lim), 'sd_BP_high_dps.mat'], 'BP_high_dps', 'BP_high_cum_dps')

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