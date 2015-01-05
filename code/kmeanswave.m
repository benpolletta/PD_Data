function Peak_data= kmeanswave
% most of the below copy pasted from PD_spikes
load('st_stn_subjects.mat')

load('13204_all_channel_data_dec.mat')

% specify the folder to use
fo = 1;
% specify the minimum prominence
min_prominences = 150;
% specify minimum separation in time between the spikes
min_secs_apart = 0.5;
% below code was just copied for simplicity
min_prominences = ones(length(folders), 1)*min_prominences;
min_secs_apart = ones(length(folders), 1)*min_secs_apart;
basetime = basetimes(fo);
% initialize time array
t = (1:length(PD_dec))/sampling_freq - basetime;
% again, the next line was just copied and pasted
min_prominence = min_prominences(fo);

min_samples_apart = min_secs_apart(fo)*sampling_freq;

spike_start = 1; spike_end = length(PD_dec);
% get LFP data from both channels
Spike_data = PD_dec(spike_start:spike_end, :);
    
Peak_data = cell(1, 2);

fig_name = sprintf('%s_%dto%d_%.1fHz_%.1fprom_spikes', '130204', round(floor(spike_start/(sampling_freq*60))), round(ceil(spike_end/(sampling_freq*60))),...
        sampling_freq/min_samples_apart, min_prominence);

    for ch = 1:2
        
        close all
        
%         find the peaks and their locations within the data
        [peaks, locs] = findpeaks(Spike_data(:, ch), 'MinPeakDistance', min_samples_apart);
        
%         get peak_details. This function was modified. 
        [peak_widths, peak_prominences, peak_forms] = peak_details(Spike_data(:, ch), locs, min_samples_apart, sampling_freq);
        
        locs = locs + spike_start - 1;
        
        peak_prom_std = zscore(peak_prominences);
        
        opts = statset('Display','final');
        %         meat of the code. Replicates ten times (picks a random
        %         start point). 
        [IDX, C] = kmeans(peak_forms,3,'Replicates',10,'Options',opts);
        
        Peak_data{ch} = [peaks locs peak_widths peak_prom_std peak_prominences IDX];

    end
   
   plot_peaks_2(fig_name, t, sampling_freq, PD_dec, Peak_data, [spike_start spike_end] + [-min_samples_apart min_samples_apart])

end

% copy pasta from PD_spikes. edited (see comments)
function [peak_widths, peak_prominences, peak_forms] = peak_details(data, locs, min_distance, sf)
% specify sampling frequency
sampling_freq = sf;

data_length = length(data);
% number of peaks is equal to the number of peak locations
no_peaks = length(locs);
% initialize peak width and peak prominence arrays
[peak_widths, peak_prominences] = deal(nan(size(locs)));
% peak forms is new. This is to store the amplitudes for the samples
% surrounding each peak
[peak_forms] = nan(no_peaks,sampling_freq);

for p = 1:no_peaks
%     get location of current peak
    loc = locs(p);
%   get beginning and end of window that surrounds the peak
    window_start = max(loc - ceil(min_distance/2), 1);
    
    window_end = min(loc + ceil(min_distance/2), data_length);
    
    window_indices = window_start:window_end;
    
%     get the difference between the peak height and the heights of the
%     surrounding sample
    dropoff = data(loc) - data(window_indices);
    
%     this will mark the beginning of the "peak" period
    peak_start = max(floor(loc-sampling_freq/2),1);

%     this will mark the end of the "peak" period
    peak_end = min(floor(loc+sampling_freq/2),length(data));
    
%     this will store the vector 
    peak_vector = data(peak_start:(peak_end-1));
    
%     subtract out mean
    peak_vector = peak_vector - mean(peak_vector);
    
%     add nans if the peak is very close to the beginning or end
    if (loc-sampling_freq/2 < 1)
        peak_vector = [nan(abs(loc-sampling_freq/2),1); peak_vector];
    end
    
    if (length(peak_vector) < sampling_freq)
        peak_vector = [peak_vector; nan(sampling_freq-length(peak_vector),1)];
    end
    
    peak_forms(p,:) = peak_vector';
    
    peak_prominences(p) = sum(abs(dropoff));
    
    d_dropoff = diff(dropoff);

    left_endpoint = find(sign(d_dropoff(window_indices < loc)) < 0, 1, 'last');
    
    offset = length(window_indices)-1-length(d_dropoff(window_indices(1:(end - 1)) > loc));
    
    right_endpoint = find(sign(d_dropoff(window_indices(1:(end - 1)) > loc)) > 0, 1)+offset;
    
    try
        peak_widths(p) = right_endpoint - left_endpoint + 1;
        assert(peak_widths(p) >= 0, ['Peak width: ',peak_widths(p)]);
    catch err
        peak_widths(p) = NaN;
        err.stack
    end
    
end

end

function plot_peaks_2(fig_name, t, sampling_freq, PD_dec, Peak_data, plot_bounds)

epoch_secs = 20;

epoch_length = epoch_secs*sampling_freq;

colors = [0 0 0; 1 0 0; 0 0 1; 0 0.75 0];

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
        
        epoch_prominences{ch} = Peak_data{ch}(epoch_loc_indices, 5);
        
        epoch_prom_std{ch} = Peak_data{ch}(epoch_loc_indices, 4);
        
        plot_peaks(epoch_locs, ch) = ((-1)^(ch + 1))*epoch_peaks{ch};
        cluster{ch} = Peak_data{ch}(epoch_loc_indices, 6);
    end
    
    if mod(e, 3) == 1

        figure

    end

    subplot(3, 1, mod(e, 3) + 3*(mod(e, 3) == 0))
    
    plot(epoch_t, epoch_data)
    
    hold on
    
    plot(epoch_t, plot_peaks, 'v')

    for ch = 1:2
%     for ch = 2
        uniques = unique(cluster{ch});
        uniques(isnan(uniques)) = [];
        for ind = 1:length(uniques)
            clust = uniques(ind);
            curr_indices = find(cluster{ch} == clust);
%             text(epoch_loc_t{ch}(curr_indices), ((-1)^(ch + 1))*epoch_peaks{ch}(curr_indices), cellstr(num2str(round(10*epoch_prominences{ch}(curr_indices))/10)),...
%                 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'Color', curr_colors(mod(clust,size(curr_colors,1))+1, :), 'FontSize', 12, 'FontWeight', 'bold')

            text(epoch_loc_t{ch}(curr_indices), ((-1)^(ch + 1))*epoch_peaks{ch}(curr_indices)+mod(ch,2)*.2, cellstr(num2str(round(10*clust)/10)),...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', 'Color', colors(clust, :), 'FontSize', 12, 'FontWeight', 'bold')
        end
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