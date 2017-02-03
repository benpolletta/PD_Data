function increase_summary = power_post_vs_pre_summary(subject_mat_cell, save_name, channel_labels, peak_suffix, norm, freqs, no_cycles, bands, band_indices, window_length)

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

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(channel_labels);

for s = 1:length(subject_mat_cell)
    
    load(subject_mat_cell{s})
    
    subject_mat_folders(s) = length(folders);
    
end

no_folders = sum(subject_mat_folders); % - 2;

no_blocks = 3;

stat_labels = {'Total Duration', 'Start Time', 'End Time'}; no_stats = length(stat_labels);

increase_summary = nan(no_stats, no_folders, no_blocks);

blocks = cell(no_blocks, 1);

folder_index = 0;

for s = 1:length(subject_mat_cell)
    
    load(subject_mat_cell{s})
    
    chan_order = [1 2];
    
    if ~strcmp(chan_labels{1}, 'Striatum');
        
        chan_order = fliplr(chan_order);
        
    end
    
    for fo = 1:subject_mat_folders(s)
        
        folder_index = folder_index + 1;
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        if ~(strcmp(folder, '130716') || strcmp(folder, '130830'))
            
            for b = band_indices
                
                load([subj_name, BP_suffix, win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_stats.mat'])
                
                for ch = 1:no_chans
                    
                    increase_blocks = index_to_blocks(test_win(:, 2, chan_order(ch)));
                    
                    increase_sec_blocks = win_secs(increase_blocks);
                    
                    blocks{ch} = increase_sec_blocks;
                    
                end
                
                overlap_index = test_win(:, 2, 1) & test_win(:, 2, 2);
                
                overlap_blocks = index_to_blocks(overlap_index);
                
                overlap_sec_blocks = win_secs(overlap_blocks);
                
                blocks{3} = overlap_sec_blocks;
                
                for block = 1:no_blocks
                    
                    if ~isempty(blocks{block})
                        
                        if size(blocks{block}, 2) ~= 2
                            
                            blocks{block} = blocks{block}';
                            
                        end
                        
                        increase_sec_blocks = blocks{block} + (window_length/2)*ones(size(blocks{block}))*diag([-1 1]);
                        
                        increase_index_secs = blocks_to_index(increase_sec_blocks, t_sec{2});
                        
                        increase_sec_blocks = index_to_blocks(increase_index_secs);
                        
                        increase_length = sum(diff(increase_sec_blocks, [], 2)); % + epoch_secs);
                        
                        increase_start = increase_sec_blocks(1, 1); % - epoch_secs/2;
                        
                        increase_end = increase_sec_blocks(end, end); % + epoch_secs/2;
                        
                        increase_summary(:, folder_index, block) = [increase_length; increase_start; increase_end]/60;
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

increase_summary = permute(increase_summary, [3 2 1]);

save([save_name, BP_suffix, norm, win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_summary.mat'], 'increase_summary')

figure

for stat = 1:no_stats
    
    h = subplot(2, no_stats, stat);
        
    set(h, 'NextPlot', 'add', 'ColorOrder', distinguishable_colors(no_folders), 'FontSize', 14)
    
    non_nan_indices = any(~isnan(increase_summary(:, :, stat)));

    plot((1:no_blocks)', increase_summary(:, non_nan_indices, stat), 'LineWidth', 2)
    
    hold on
    
    for block = 1:no_blocks
    
        plot(block, increase_summary(block, non_nan_indices, stat)', 'o', 'MarkerSize', 8, 'LineWidth', 1)
        
    end
    
    set(gca, 'XLim', [.5 no_stats + .5], 'XTick', 1:no_stats, 'XTickLabel', {channel_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
    title(stat_labels{stat})
    
    subplot(2, no_stats, no_stats + stat)
    
    h = barwitherr(nanstd(increase_summary(:, :, stat)'), nanmean(increase_summary(:, :, stat)'));
    
    set(h, 'FaceColor', [.5 .5 .5])
    
    set(gca, 'XLim', [.25 no_stats + .75], 'XTickLabel', {channel_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
end

save_as_pdf(gcf, [save_name, BP_suffix, norm, win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_summary'])

end
