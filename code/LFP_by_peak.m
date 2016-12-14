function LFP_by_peak(subj_name, chan_labels)

if isempty(chan_labels)
    
    chan_labels = cell(2, 1);
    
    for ch = 1:2
        
        chan_labels{ch} = sprintf('Ch. %d', ch);
        
    end
    
end

pair_colors = {[0 0 0; 0 0 1], [0 0 0; 0 .5 0]};

load([subj_name, '_all_channel_data_dec.mat'])

[Spike_indicator, ~] = peak_loader(subj_name, '_kmeans', length(PD_dec));

chan_pairs = [2 1; 1 2];

chan_linewidth = [1 2];

figure

%% Subtracting peaks.

for p = 1:2
    
    clear Spike_locs
    
    Spike_locs = find(Spike_indicator(:, chan_pairs(p, 1)));
    
    min_peak_distance = max(min(diff(Spike_locs)), sampling_freq);
    
    for ch = 1:2
        
        clear data sum_peak_form
        
        data = PD_dec(:, chan_pairs(p, ch));
        
        peak_forms = get_peak_forms(data, Spike_locs, min_peak_distance);
        
        sum_peak_form = nansum(peak_forms);
        
        time_re_peak = ((1:length(sum_peak_form)) - round(length(sum_peak_form)/2))/sampling_freq;
        
        subplot(1,2,p)
        
        plot(time_re_peak, sum_peak_form, 'Color', pair_colors{p}(ch, :), 'LineWidth', chan_linewidth(ch))
        
        axis tight, hold on
        
        legends{ch} = sprintf('%s LFP by %s Peaks', chan_labels{chan_pairs(p, ch)}, chan_labels{chan_pairs(p, 1)});
        
    end
    
    legend(legends, 'Location', 'Southwest')
    
    y_limits = ylim;
    
    plot([0 0]', y_limits', 'k')
    
    title(sprintf('%s, %s LFP (Mean) Relative to %s Peaks', subj_name, chan_labels{fliplr(chan_pairs(p, :))}), 'FontSize', 16)
    
    xlabel(sprintf('Time Rel. %s Peak Location (s)', chan_labels{chan_pairs(p, 1)}), 'FontSize', 16)
    
    ylabel(sprintf('Average %s LFP', chan_labels{chan_pairs(p, 2)}), 'FontSize', 16)
        
end

save_as_pdf(gcf, [subj_name, '_LFP_by_peak'])

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