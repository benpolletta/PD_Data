function data_subtracted = subtract_peaks(subj_name, peak_suffix, outlier_lim, plot_opt)

load([subj_name, '_all_channel_data_dec.mat'])

data_subtracted = PD_dec;

[Spike_indicator, Peak_data] = peak_loader(subj_name, peak_suffix, length(PD_dec));

data_subtracted = PD_dec;

for ch = 1:2
    
    data = PD_dec(:, ch);
    
    if isempty(Peak_data{ch}) | isnan(Peak_data{ch})
        
        Spike_locs = find(Spike_indicator(:, ch));
        
        min_peak_distance = min(diff(Spike_locs));
        
        peak_forms = get_peak_forms(data, Spike_locs, min_peak_distance);
        
        opts = statset('Display','final');
        
        no_clusters = 5; clusters = 1:no_clusters;
        
        [IDX, ~] = kmeans(peak_forms, no_clusters, 'Replicates', 10, 'Options', opts);
        
        Peak_data{ch} = [nan(size(Spike_locs)) Spike_locs];
        
        Peak_data{ch}(:, 6) = IDX;
        
    else

        min_peak_distance = min(diff(Peak_data{ch}(:, 2)));
        
        [clusters, no_spikes] = find_clusters(Spike_indicator(:, ch), Peak_data{ch});
        
        no_clusters = length(clusters);
        
        if no_clusters < 5 & no_spikes > 100
            
            Peak_data{ch} = [];
            
            Spike_locs = find(Spike_indicator(:, ch));
            
            peak_forms = get_peak_forms(data, Spike_locs, min_peak_distance);
            
            opts = statset('Display','final');
            
            no_clusters = 5; clusters = 1:no_clusters;
            
            [IDX, ~] = kmeans(peak_forms, no_clusters, 'Replicates', 10, 'Options', opts);
        
            Peak_data{ch} = [nan(size(Spike_locs)) Spike_locs];
            
            Peak_data{ch}(:, 6) = IDX;
            
        end
        
    end
    
    for c = 1:no_clusters
            
        w = 0;
        
        cluster_peak_locs = Peak_data{ch}(Peak_data{ch}(:, 6) == clusters(c), 2);
        
        no_peaks = size(cluster_peak_locs, 1);
        
        standard_window = -ceil(min_peak_distance/2):ceil(min_peak_distance/2);
        
        x_envelope = standard_window/max(standard_window);
        
        envelope = exp(-1./(1 - x_envelope.^2))';
        
        envelope = envelope/max(envelope);
            
        peak_forms = get_peak_forms(data, cluster_peak_locs, min_peak_distance);
        
        if no_peaks > 1
            
            cluster_centroid = nanmean(peak_forms)';
            
            cluster_centroid = cluster_centroid - nanmean(cluster_centroid);
            
            cluster_centroid = cluster_centroid/norm(cluster_centroid, 2);
            
            for p = 1:no_peaks
                
                % Get location of current peak.
                loc = cluster_peak_locs(p);
                
                % Get beginning and end of window that surrounds the peak.
                window_start = max(loc - ceil(min_peak_distance/2), 1);
                window_end = min(loc + ceil(min_peak_distance/2), length(data));
                window_indices = window_start:window_end;
                window_indicator = standard_window >= (window_start - loc) & standard_window <= (window_end - loc);
                
                peak_data = data(window_indices);
                
                peak_proj = (peak_data'*cluster_centroid(window_indicator))*cluster_centroid(window_indicator);
                
                peak_subtracted = peak_data - envelope(window_indicator).*peak_proj;
                
                data_subtracted(window_indices, ch) = peak_subtracted;
                
                if plot_opt && ~w
                    
                    figure
                    
                    subplot(211) % (311)
                    
                    plot(window_indices, [peak_data, cluster_centroid(window_indicator),...
                        peak_subtracted, peak_proj, envelope(window_indicator).*peak_proj])
                    
                    axis tight
                    
                    legend({'Data', 'Centroid', 'Subtracted', 'Projection', 'Windowed Proj.'})
                    
                    title(sprintf('%s, Channel %d, Cluster %d, Peak %d', subj_name, ch, c, p))
                    
                    peak_wav = wavelet_spectrogram(peak_data, sampling_freq, 1:200, linspace(3,21,200), 0, '');
                    
                    % subplot(312)
                    %
                    % imagesc(window_indices, 1:200, nanzscore(abs(peak_wav))')
                    %
                    % axis xy
                    
                    subtracted_wav = wavelet_spectrogram(peak_subtracted, sampling_freq, 1:200, linspace(3,21,200), 0, '');
                    
                    subplot(212) % (313)
                    
                    imagesc(window_indices, 1:200, nanzscore(abs(peak_wav) - abs(subtracted_wav))') % (nanzscore(abs(peak_wav)) - nanzscore(abs(subtracted_wav)))')
                    
                    axis xy
                    
                    save_as_pdf(gcf, sprintf('%s_chan%d_clust%d_peak%d', subj_name, ch, c, p))
                    
                    w = waitforbuttonpress;
                
                    if w, close('all'), end
                    
                end
                
            end
            
        else
            
            % Get location of current peak.
            loc = cluster_peak_locs(1);
            
            % Get beginning and end of window that surrounds the peak.
            window_start = max(loc - ceil(min_peak_distance/2), 1);
            window_end = min(loc + ceil(min_peak_distance/2), length(data));
            window_indices = window_start:window_end;
            window_indicator = standard_window >= (window_start - loc) & standard_window <= (window_end - loc);
            
            peak_data = data(window_indices);
            
            trend = polyfit(window_indices', peak_data, 1);
            
            peak_trend = trend(1)*window_indices' + trend(2);
            
            peak_detrended = peak_data - peak_trend;
            
            peak_subtracted = peak_data - envelope(window_indicator).*peak_detrended;
            
            data_subtracted(window_indices, ch) = peak_subtracted;
            
            if plot_opt && ~w
                
                figure
                
                subplot(211) % (311)
                
                plot(window_indices, [peak_data, peak_detrended,...
                    peak_subtracted, envelope(window_indicator).*peak_detrended])
                
                axis tight
                
                legend({'Data', 'Detrended', 'Subtracted', 'Windowed'})
                
                title(sprintf('%s, Channel %d, Cluster %d, Peak 1', subj_name, ch, c))
                
                peak_wav = wavelet_spectrogram(peak_data, sampling_freq, 1:200, linspace(3,21,200), 0, '');
                
                % subplot(312)
                %
                % imagesc(window_indices, 1:200, nanzscore(abs(peak_wav))')
                %
                % axis xy
                
                subtracted_wav = wavelet_spectrogram(peak_subtracted, sampling_freq, 1:200, linspace(3,21,200), 0, '');
                
                subplot(212) % (313)
                
                imagesc(window_indices, 1:200, nanzscore(abs(peak_wav) - abs(subtracted_wav))') % (nanzscore(abs(peak_wav)) - nanzscore(abs(subtracted_wav)))')
                
                axis xy
                
                save_as_pdf(gcf, sprintf('%s_chan%d_clust%d_peak1', subj_name, ch, c))
                
                w = waitforbuttonpress;
                
                if w, close('all'), end
                
            end
            
        end
        
    end
    
end

save(sprintf('%s_all_channel_data_dec_peakless', subj_name), 'data_subtracted')

end

function [clusters, total_spikes] = find_clusters(Spike_indicator, Peak_data)

clusters = []; total_spikes = 0;

cluster_ids = Peak_data(:, 6);

cluster_list = unique(cluster_ids);

cluster_list(isnan(cluster_list)) = [];

for i = cluster_list'
    
    cluster_indicator = zeros(size(Spike_indicator));
    
    cluster_indicator(Peak_data(cluster_ids == i, 2)) = 1;
    
    non_cluster_spikes = max(Spike_indicator - cluster_indicator, 0);
    
    if sum(non_cluster_spikes) < sum(Spike_indicator)
       
        clusters = [clusters i];
        
        total_spikes = total_spikes + sum(cluster_indicator);
        
    end
    
end

end

function peak_forms = get_peak_forms(data, locs, min_peak_distance)
    
    data_length = length(data);
    
    no_peaks = length(locs);
    
    standard_window = -ceil(min_peak_distance/2):ceil(min_peak_distance/2);
    
    peak_forms = nan(no_peaks, length(standard_window));

    for p = 1:no_peaks

        % Get location of current peak.
        loc = locs(p);
        
        % Get beginning and end of window that surrounds the peak.
        window_start = max(loc - ceil(min_peak_distance/2), 1);
        window_end = min(loc + ceil(min_peak_distance/2), data_length);
        window_indices = window_start:window_end;
        window_indicator = standard_window >= (window_start - loc) & standard_window <= (window_end - loc);

        % Get peak.
        peak_forms(p, window_indicator) = data(window_indices);

    end

end



% %%
% 
% if any(Spike_indicator ~= 0)
%     
%     [spike_nans_wt, spike_nans] = indicator_to_nans(Spike_indicator, sampling_freq, freqs, no_cycles, bands);
%     
%     BP(logical(spike_nans)) = nan;
%     
%     Spec(logical(spike_nans_wt)) = nan;
%     
% end
% 
% if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
%     
%     load([subj_name, '_wav_laser_artifacts.mat'])
%     
%     [laser_nans_wt, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, no_cycles, bands);
%     
%     laser_nans = repmat(laser_nans, [1 1 2]);
%     
%     BP(logical(laser_nans)) = nan;
%     
%     laser_nans_wt = repmat(laser_nans_wt, [1 1 2]);
%     
%     Spec(logical(laser_nans_wt)) = nan;
%     
% end
% 
% if ~isempty(outlier_lim) && ~isempty(dir([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat']))
%     
%     load([subj_name, '_wav_BP_', num2str(outlier_lim), 'sd_outliers.mat'])
%     
%     [outlier_nans_wt, outlier_nans] = indicator_to_nans(double(artifact_indicator), sampling_freq, freqs, no_cycles, bands);
%     
%     outlier_nans = repmat(outlier_nans, [1 1 2]);
%     
%     BP(logical(outlier_nans)) = nan;
%     
%     outlier_nans_wt = repmat(outlier_nans_wt, [1 1 2]);
%     
%     Spec(logical(outlier_nans_wt)) = nan;
%     
% end