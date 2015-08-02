function PD_spikes_by_freq_individual_save(folder, prefix, channel_multipliers, spike_boundaries, band, no_cycles, min_secs_apart, min_prominence, min_prominences_flag)

% Script to save peaks in time series from carbachol data. To help in
% locating peaks, this script first finds peaks in band power, then locates
% peaks in the LFP (above a certain prominence cutoff) that are close to
% these peaks in band power.

% INPUTS:
% folder - folder name (string).
% prefix - prefix name (string).
% basetime - time of infusion/first laser on, in seconds (integer).
% channel_multipliers - searches for peaks in LFP if +1, troughs in LFP if
%   -1 (matrix, length = number of channels).
% spike_boundaries - times within which to search for peaks/troughs
%   (in samples, integers).
% band - limits (in Hz) of band close to whose peaks spikes are sought, e.g. [1 4] (length 2 vector).
% band_label - name of band for saving figures, e.g. 'delta' (string).
% no_cycles - number of cycles to use for Morlet wavelets used to get band
%   power, linearly spaced from lowest frequency up to highest frequency, e.g. [3 4] (length 2 vector). 
% min_secs_aparts - minimum number of seconds separating consecutive peaks
%   (seconds, fload).
% min_prominence - prominence cutoff at which peaks are recognized as
% peaks. Prominence is the sum over a certain time window (given by
%   min_secs_apart) of the difference between the peak value and the
%   surrounding values.
% min_prominence_flag - 'std', if prominence cutoff is given in terms of 
%   z-scored prominence for all peaks (determined before peaks are excluded
%   by prominence), or '' if prominence cutoff is given in terms of raw
%   prominence (usually easier).

subj_name = [folder, '/', prefix];

load([subj_name, '_all_channel_data_dec.mat'])

min_samples_apart = min_secs_apart*sampling_freq;

if isempty(channel_multipliers)
   
    channel_multipliers = [1 1];
    
end

if isempty(spike_boundaries)
    
    spike_start = 1; spike_end = length(PD_dec);
    
else
    
    spike_start = spike_boundaries(1);
    
    spike_end = spike_boundaries(2);
    
end

% fig_name = sprintf('%s_%dto%d_%.1fHz_%.1f%sprom_spikes', subj_name, round(floor(spike_start/(sampling_freq*60))), round(ceil(spike_end/(sampling_freq*60))),...
%     sampling_freq/min_samples_apart, min_prominence, min_prominences_flag);

%% Finding spikes, their width and prominence.

Spike_data = PD_dec(spike_start:spike_end, :);

Spike_data_wav = wavelet_spectrogram(Spike_data, sampling_freq, band(1):band(2), no_cycles(1):no_cycles(2), 0, '');

Spike_data_BP = permute(sum(abs(Spike_data_wav), 2), [1 3 2]);

Spike_indicator = zeros(size(PD_dec));

Peak_data = cell(1, 2);

for ch = 1:2
    
    [~, BP_locs] = findpeaks(Spike_data_BP(:, ch), 'MinPeakDistance', min_samples_apart); %((-1)^(ch + 1))*
    
    [peaks, locs] = BP_to_LFP(channel_multipliers(ch)*Spike_data(:, ch), BP_locs, 100);
    
    [peak_widths, peak_prominences] = peak_details(channel_multipliers(ch)*Spike_data(:, ch), locs, min_samples_apart);
    
    locs = locs + spike_start - 1;
    
    peak_prom_std = zscore(peak_prominences);
    
    Peak_data{ch} = [peaks locs peak_widths peak_prom_std peak_prominences];
    
    if strcmp(min_prominences_flag, '')
        
        Peak_data{ch}(peak_prominences < min_prominence(ch), :) = [];
        
    elseif strcmp(min_prominences_flag, 'std')
        
        Peak_data{ch}(peak_prom_std < min_prominence(ch), :) = [];
        
    end
    
    Spike_indicator(Peak_data{ch}(:, 2), ch) = 1;
    
end

%% Saving peaks.

save([subj_name, '_peaks.mat'], 'Peak_data', 'Spike_indicator', 'band', 'no_cycles', 'min_secs_apart', 'min_prominence', 'min_prominences_flag')

end

function [LFP_peaks, LFP_peak_locs] = BP_to_LFP(LFP, BP_peaks, window)

no_peaks = length(BP_peaks);

[LFP_peaks, LFP_peak_locs] = deal(nan(size(BP_peaks)));

LFP_length = length(LFP);

for p = 1:no_peaks
   
    BP_peak_loc = BP_peaks(p);
    
    LFP_start = max(round(BP_peak_loc - window/2), 1);
    
    LFP_end = min(round(BP_peak_loc + window/2), LFP_length);
    
    LFP_around_peak = LFP(LFP_start:LFP_end);
    
    [peak, loc] = findpeaks(LFP_around_peak, 'MinPeakDistance', length(LFP_around_peak) - 1);
    
    LFP_peaks(p) = peak;
    
    LFP_peak_locs(p) = loc + LFP_start - 1;
    
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

