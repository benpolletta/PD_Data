function WT_by_peak(subj_name, chan_labels)

if isempty(chan_labels)
    
    chan_labels = cell(2, 1);
    
    for ch = 1:2
        
        chan_labels{ch} = sprintf('Ch. %d', ch);
        
    end
    
end

pair_colors = {[0 0 0; 0 0 1], [0 0 0; 0 .5 0]};

load([subj_name, '_all_channel_data_dec.mat'])

[Spike_indicator, ~] = peak_loader(subj_name, '_kmeans', length(PD_dec));

chan_pairs = [1 2; 2 1];

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
        
        [avg_wt_power, ~, time_re_peak] = get_peak_triggered_wt(peak_forms, 1:100, linspace(3, 10, 100), sampling_freq);

        peak_loc = zeros(size(time_re_peak))'; peak_loc(time_re_peak == min(abs(time_re_peak))) = 1;
        
        [Spec_nans, ~] = indicator_to_nans(peak_loc, sampling_freq, 1:100, linspace(3, 10, 100), []);
        
        subplot(2, 2, (ch - 1)*2 + p)
        
        avg_wt_power(logical(Spec_nans)') = nan;
        
        imagesc(time_re_peak, 1:100, nanzscore(avg_wt_power')')
        
        axis xy, hold on
        
        colorbar
        
        plot([0 0]', [1 100]', 'w', 'LineWidth', 1.5)
        
        plot(time_re_peak([1 end])', [15 30; 15 30], 'w', 'LineWidth', 1.5)
        
        title({sprintf('%s, %s Wavelet Transform (Mean)', subj_name, chan_labels{chan_pairs(p, ch)});...
            sprintf('Relative to %s Peaks', chan_labels{chan_pairs(p, 1)})}, 'FontSize', 16)
        
        xlabel(sprintf('Time Rel. %s Peak Location (s)', chan_labels{chan_pairs(p, 1)}), 'FontSize', 16)
        
        ylabel('Freq. (Hz)', 'FontSize', 16)
        
    end
        
end

save_as_pdf(gcf, [subj_name, '_WT_by_peak'])

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

function [avg_wt_power, avg_wt_phase, time_re_peak] = get_peak_triggered_wt(peak_forms, freqs, no_cycles, sampling_freq)

data_size = size(peak_forms, 2);

no_freqs = length(freqs);

cycle_lengths = no_cycles.*(sampling_freq./freqs);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', max(sampling_freq, max(cycle_lengths)));

peak_loc = zeros(data_size, 1); peak_loc(round(data_size/2)) = 1;

[Spec_nans, ~] = indicator_to_nans(peak_loc, sampling_freq, freqs, no_cycles, []); Spec_nans = Spec_nans';

segment_length = sampling_freq;

avg_wt_power = nan(no_freqs, data_size);

avg_wt_phase = nan(no_freqs, data_size);

no_peaks = size(peak_forms, 1);

wt_power = nan(no_freqs, data_size, no_peaks);

wt_phase = nan(no_freqs, data_size, no_peaks);

parfor p = 1:no_peaks
    
    Spec = nan(no_freqs, data_size);
    
    peak_form = peak_forms(p, :);
    
    peak_form_reflected = [fliplr(peak_form(1:segment_length)), peak_form, fliplr(peak_form((end - segment_length + 1):end))];
    
    for f = 1:no_freqs
        
        conv_prod = nanconv(peak_form_reflected, wavelets(f,:), 'same');
        
        Spec(f, :) = conv_prod((segment_length + 1):(end - segment_length));
        
    end
    
    % Spec(logical(Spec_nans)) = nan;
    
    wt_power(:, :, p) = abs(Spec);
    
    wt_phase(:, :, p) = angle(Spec);
    
end

avg_wt_power(:, :) = nanmean(wt_power, 3);

avg_wt_phase(:, :) = nanmean(wt_phase, 3);

time_re_peak = ((1:data_size) - data_size/2)/sampling_freq;

end