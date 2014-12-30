function PD_spikes(subjects_mat, spike_boundaries, min_secs_apart, min_prominences, min_prominences_flag, plot_opt)

% Script to find spikes in time series from carbachol data.

load(subjects_mat)

if length(min_secs_apart) == 1
    
    min_secs_apart = ones(length(folders), 1)*min_secs_apart;
    
end

if length(min_prominences) == 1
    
    min_prominences = ones(length(folders), 1)*min_prominences;
    
end

load([folders{1}, '/', prefixes{1}, '_all_channel_data_dec.mat'], 'sampling_freq')

for fo = 1:1%length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    basetime = basetimes(fo);
    
    t = (1:length(PD_dec))/sampling_freq - basetime;
    
    min_prominence = min_prominences(fo);
    
    min_samples_apart = min_secs_apart(fo)*sampling_freq;
    
    if isempty(spike_boundaries)
        
        spike_start = 1; spike_end = length(PD_dec);
        
    else
        
        spike_start = spike_boundaries(fo, 1);
        
        spike_end = spike_boundaries(fo, 2);
        
    end
    
    fig_name = sprintf('%s_%dto%d_%.1fHz_%.1f%sprom_spikes', subj_name, round(floor(spike_start/(sampling_freq*60))), round(ceil(spike_end/(sampling_freq*60))),...
        sampling_freq/min_samples_apart, min_prominence, min_prominences_flag);
    
    %% Finding spikes, their width and prominence.
    
    Spike_data = PD_dec(spike_start:spike_end, :);
    
    Peak_data = cell(1, 2);
    
    for ch = 1:2
        
        [peaks, locs] = findpeaks(Spike_data(:, ch), 'MinPeakDistance', min_samples_apart); %((-1)^(ch + 1))*
        
        [peak_widths, peak_prominences] = peak_details(Spike_data(:, ch), locs, min_samples_apart);
        
        % peaks(peak_prominences < min_prominence) = [];
        % 
        % locs(peak_prominences < min_prominence) = [];
        % 
        % peak_prominences(peak_prominences < min_prominence) = [];
        % 
        % peak_widths(peak_prominences < min_prominence) = [];
        
        locs = locs + spike_start - 1;
        
        peak_prom_std = zscore(peak_prominences);
        
        Peak_data{ch} = [peaks locs peak_widths peak_prom_std peak_prominences];
        
        if strcmp(min_prominences_flag, '')
        
            Peak_data{ch}(peak_prominences < min_prominence, :) = [];
            
        elseif strcmp(min_prominences_flag, 'std')
        
            Peak_data{ch}(peak_prom_std < min_prominence, :) = [];
            
        end
            
    end
    
    %% Plotting peaks.
    
    if plot_opt == 1
        
        plot_peaks_1(t, PD_dec, Peak_data, min_samples_apart)
    
    elseif plot_opt == 2
    
        plot_peaks_2(fig_name, t, sampling_freq, PD_dec, Peak_data, [spike_start spike_end] + [-min_samples_apart min_samples_apart])
        
    end
        
end

end

function [peak_widths, peak_prominences] = peak_details(data, locs, min_distance)

data_length = length(data);

no_peaks = length(locs);

[peak_widths, peak_prominences] = deal(nan(size(locs)));

for p = 1:no_peaks
    
    loc = locs(p);
    
    window_start = max(loc - ceil(min_distance/2), 1);
    
    window_end = min(loc + ceil(min_distance/2), data_length);
    
    window_indices = window_start:window_end;
    
    dropoff = data(loc) - data(window_indices);
    
    peak_prominences(p) = sum(abs(dropoff));
    
    d_dropoff = diff(dropoff);
    
    left_endpoint = find(sign(d_dropoff(window_indices < loc)) < 0, 1, 'last');
    
    right_endpoint = find(sign(d_dropoff(window_indices(1:(end - 1)) > loc)) > 0, 1);
    
    peak_widths(p) = right_endpoint - left_endpoint + 1;
    
end

end

function plot_peaks_1(t, PD_dec, Peak_data, min_sample_apart)

no_rows = 8;

colors = {'b', 'g'};

symbols = {'v', '^'};

for ch = 1:2
    
    peaks = Peak_data{ch}(:, 1); locs = Peak_data{ch}(:, 2); peak_prominences = Peak_data{ch}(:, end);
    
    no_peaks = length(peaks);
    
    no_figs = floor(no_peaks/(no_rows^2));
    
    p = 1;
    
    for f = 1:no_figs
        
        figure
        
        for s = 1:(no_rows^2)
            
            subplot(no_rows, no_rows, s)
            
            peak_loc = locs(p);
            
            plot_span = (peak_loc - 3*min_sample_apart):(peak_loc + 3*min_sample_apart);
            
            plot(t(plot_span), PD_dec(plot_span, :))
            
            hold on
            
            plot(t(locs(p)), ((-1)^(ch + 1))*peaks(p), [colors{ch}, symbols{ch}])
            
            title(num2str(peak_prominences(p)))
            
            axis tight
            
            p = p + 1;
            
        end
        
    end
    
    figure
    
    no_leftovers = no_peaks - no_figs*(no_rows^2);
    
    [r, c] = subplot_size(no_leftovers);
    
    l = 1;
    
    for s = 1:no_leftovers
        
        subplot(r, c, s)
        
        peak_loc = locs(no_figs*(no_rows^2) + l);
        
        plot_span = (peak_loc - 3*min_sample_apart):(peak_loc + 3*min_sample_apart);
        
        plot(t(plot_span), PD_dec(plot_span, :))
        
        hold on
        
        plot(t(locs(no_figs*(no_rows^2) + l)), ((-1)^(ch + 1))*peaks(no_figs*(no_rows^2) + l), [colors{ch}, symbols{ch}])
        
        title(num2str(peak_prominences(no_figs*(no_rows^2) + l)))
        
        axis tight
        
        l = l + 1;
        
    end
    
end

end

function plot_peaks_2(fig_name, t, sampling_freq, PD_dec, Peak_data, plot_bounds)

epoch_secs = 20;

epoch_length = epoch_secs*sampling_freq;

colors = [.75 0 0; 1 0 0];

no_epochs = ceil((diff(plot_bounds) + 1)/epoch_length);

for e = 1:no_epochs

    epoch_start = max(plot_bounds(1) + (e - 1)*epoch_length, 1);

    epoch_end = min(plot_bounds(1) + e*epoch_length - 1, length(PD_dec));
    
    epoch_t = t(epoch_start:epoch_end);

    epoch_data = PD_dec(epoch_start:epoch_end, :);
    
    plot_peaks = nan(length(epoch_t), 2);
    
    [epoch_loc_t, epoch_peaks, epoch_prominences] = deal(cell(1, 2));
    
    for ch = 1:2
        
        epoch_loc_indices = Peak_data{ch}(:, 2) >= epoch_start & Peak_data{ch}(:, 2) <= epoch_end;
        
        epoch_locs = Peak_data{ch}(epoch_loc_indices, 2) - epoch_start + 1;
        
        epoch_loc_t{ch} = epoch_t(epoch_locs);
        
        epoch_peaks{ch} = Peak_data{ch}(epoch_loc_indices, 1);
        
        epoch_prominences{ch} = Peak_data{ch}(epoch_loc_indices, end);
        
        epoch_prom_std{ch} = Peak_data{ch}(epoch_loc_indices, end - 1);
        
        plot_peaks(epoch_locs, ch) = epoch_peaks{ch}; %((-1)^(ch + 1))*epoch_peaks{ch};
        
    end
    
    if mod(e, 3) == 1

        figure

    end

    subplot(3, 1, mod(e, 3) + 3*(mod(e, 3) == 0))
    
    plot(epoch_t, epoch_data)
    
    hold on
    
    plot(epoch_t, plot_peaks, 'v')
    
    for ch = 1:2
        
        text(epoch_loc_t{ch}, epoch_peaks{ch}, cellstr(num2str(round(10*epoch_prominences{ch})/10)),...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', colors(ch, :), 'FontSize', 12, 'FontWeight', 'bold') %((-1)^(ch + 1))*
        
        text(epoch_loc_t{ch}, epoch_peaks{ch}, cellstr(num2str(round(10*epoch_prom_std{ch})/10)),...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', colors(ch, :), 'FontSize', 12, 'FontWeight', 'bold') %((-1)^(ch + 1))*
    
    end

    axis tight
    
    title(fig_name(1:6))
    
    if mod(e, 3) == 0 || e == no_epochs
    
        title(fig_name(1:6))
       
        minute = ceil(e/3);
        
        save_as_pdf(gcf, [fig_name, '_', num2str(minute, '%03d'), 'min'])
        
    end
    
end

end

