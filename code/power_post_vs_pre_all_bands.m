function power_post_vs_pre_all_bands(peak_suffix, norm, freqs, no_cycles, bands, window_length, significance, box_length, box_minimum)

filename = 'STR_w_M1';
channel_labels = {'Striatum', 'M1'};

if nargin < 7, significance = []; end
if nargin < 8, box_length = []; end
if nargin < 9, box_minimum = []; end

if isempty(box_length), box_length = 1; end
if isempty(box_minimum), box_minimum = 1; end
if isempty(significance), significance = .05; end

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

sampling_freq = 500;

if isempty(window_length) || window_length == 150
    
    window_length = 150; 

    win_flag = '';
    
else
    
    win_flag = ['_win', num2str(window_length)];
    
end
    
if box_length == 1 && box_minimum == 1
    
    box_flag = '';
    
else
    
    box_flag = sprintf('_%gover%d', box_minimum, box_length);
    
end

if significance == .05, sig_flag = ''; else sig_flag = sprintf('_p%g', significance); end

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(channel_labels);

for b = 1:(no_bands - 1)

    load([filename, BP_suffix, norm, win_flag, '_', short_band_labels{b}, box_flag, sig_flag, '_power_post_vs_pre_summary.mat'], 'increase_summary')
    
    all_increase_summary(:, :, :, b) = increase_summary;
    
end

stat_titles = {'Total Duration', 'Start Time', 'End Time'};
stat_labels = {'total', 'start', 'end'};
no_stats = length(stat_labels);


for stat = 1:no_stats
    
    figure
    
    barwitherr(squeeze(nanstd(all_increase_summary(:, :, stat, 1:6), [], 2)), squeeze(nanmean(all_increase_summary(:, :, stat, 1:6), 2)))
    
    set(gca, 'XLim', [.5 no_stats + .5], 'XTick', 1:no_stats, 'XTickLabel', {channel_labels{:}, 'Overlap'}, 'FontSize', 16)
    
    legend(short_band_labels(1:6))
    
    ylabel('Time (min.)')
    
    title([stat_titles{stat}, ' of Band Power Increases'])
    
    save_as_pdf(gcf, [filename, BP_suffix, norm, win_flag, '_', short_band_labels{b}, box_flag, sig_flag, '_power_post_vs_pre_all_bands_', stat_labels{stat}, '_summary'])
    
end
