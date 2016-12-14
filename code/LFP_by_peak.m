function LFP_by_peak(subj_name)

% format long g

pair_colors = [0 0 1; 0 .5 0];

load([subj_name, '_all_channel_data_dec.mat'])

[Spike_indicator, ~] = peak_loader(subj_name, '_kmeans', length(PD_dec));

chan_pairs = [2 1; 1 2];

figure

%% Subtracting peaks.

for p = 1:2
    
    data = PD_dec(:, chan_pairs(p, 2));
    
    Spike_locs = find(Spike_indicator(:, chan_pairs(p, 1)));
    
    min_peak_distance = max(min(diff(Spike_locs)), sampling_freq);
    
    peak_forms = get_peak_forms(data, Spike_locs, min_peak_distance);
    
    sum_peak_form = nansum(peak_forms);
    
    time_re_peak = ((1:length(sum_peak_form)) - round(length(sum_peak_form)/2))/sampling_freq;
    
    subplot(1,2,p)
    
    plot(time_re_peak, sum_peak_form, 'Color', pair_colors(p, :))
    
    axis tight, hold on
    
    y_limits = ylim;
    
    plot([0 0]', y_limits', 'k')
    
    title(sprintf('Channel %d LFP (Mean) Relative to Channel %d Peaks', fliplr(chan_pairs(p, :))), 'FontSize', 16)
    
    xlabel(sprintf('Time Rel. Channel %d Peak Location (s)', chan_pairs(p, 1)), 'FontSize', 16)
    
    ylabel(sprintf('Average Channel %d LFP', chan_pairs(p, 2)), 'FontSize', 16)
        
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