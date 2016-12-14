function peak_by_peak_hist(subj_name)

% format long g

pair_colors = [0 0 1; 0 .5 0];

load([subj_name, '_all_channel_data_dec.mat'])

[Spike_indicator, ~] = peak_loader(subj_name, '_kmeans', length(PD_dec));

chan_pairs = [1 2; 2 1];

figure

%% Subtracting peaks.

for p = 1:2
    
    data = Spike_indicator(:, chan_pairs(p, 2));
    
    Spike_locs = find(Spike_indicator(:, chan_pairs(p, 1)));
    
    min_peak_distance = max(min(diff(Spike_locs)), sampling_freq);
    
    peak_forms = get_peak_forms(data, Spike_locs, min_peak_distance);
    
    sum_peak_form = nansum(peak_forms); peak_length = length(sum_peak_form);
    
    time_re_peak = ((1:length(sum_peak_form)) - peak_length/2)/sampling_freq;
    
    % no_bins = 50;
    % 
    % bin_width = round(peak_length/no_bins);
    % 
    % bin_edges = ((-(bin_width/2):bin_width:(peak_length + bin_width/2)) - peak_length/2)/sampling_freq;
    % 
    % peak_count = nan(no_bins + 1, 1);
    % 
    % for b = 1:(no_bins + 1)
    % 
    %         bin_id = time_re_peak >= bin_edges(b) & time_re_peak < bin_edges(b + 1);
    % 
    %         bin_centers(b) = (bin_edges(b) + bin_edges(b + 1))/2;
    % 
    %         peak_count(b) = nansum(sum_peak_form(bin_id));
    % 
    % end
    % 
    % peak_count = peak_count/(nansum(peak_count));
    
    hold on
    
    plots(p) = plot(time_re_peak, sum_peak_form, 'Color', pair_colors(p, :));
    
    axis tight
    
    legends{p} = sprintf('Channel %d Peaks Relative to Channel %d Peaks', fliplr(chan_pairs(p, :)));
    
    % subplot(2,2, 2 + p)
    % 
    % bar(bin_centers, peak_count)
    % 
    % title(sprintf('Channel %d Peaks Relative to Channel %d Peaks', fliplr(chan_pairs(p, :))))
    % 
    % xlabel(sprintf('Time Rel. Channel %d Peak Location', chan_pairs(p, 1)))
    % 
    % ylabel('Peak Count')
        
end

legend(plots, legends)

y_limits = ylim;

plot([0 0]', y_limits', 'k:')

title('Peak-by-Peak Histograms', 'FontSize', 16)

xlabel(sprintf('Time Rel. Peak Location (s)', chan_pairs(p, 1)), 'FontSize', 16)

ylabel('Peak Count', 'FontSize', 16)

save_as_pdf(gcf, [subj_name, '_peak_by_peak_hist'])

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