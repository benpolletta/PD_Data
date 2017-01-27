function increase_summary = power_post_vs_pre_summary(subject_mat, peak_suffix, epoch_secs, pd_handle, norm, freqs, no_cycles, bands, band_indices)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

load(subject_mat)

sampling_freq = 500;

no_folders = length(folders);

if isempty(window_length)
    
    window_length = epoch_secs; 

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

no_chans = length(chan_labels);

no_blocks = 3;

stat_labels = {'Total Duration', 'Start Time', 'End Time'}; no_stats = length(stat_labels);

increase_summary = nan(no_stats, no_folders, no_blocks);

blocks = cell(3, 1);

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for b = band_indices
            
        load([subj_name, BP_suffix,win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_stats.mat'])
        
        for ch = 1:no_chans
            
            increase_blocks = index_to_blocks(test_win(:, 2, ch));
            
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
                
                increase_sec_blocks = blocks{block} + (epoch_secs/2)*ones(size(blocks{block}))*diag([-1 1]);
                
                increase_index_secs = blocks_to_index(increase_sec_blocks, t_sec{2});
                
                increase_sec_blocks = index_to_blocks(increase_index_secs);
                
                increase_length = sum(diff(increase_sec_blocks, [], 2)); % + epoch_secs);
                
                increase_start = increase_sec_blocks(1, 1); % - epoch_secs/2;
                
                increase_end = increase_sec_blocks(end, end); % + epoch_secs/2;
                
                increase_summary(:, fo, block) = [increase_length; increase_start; increase_end]/60;
                
            end
            
        end
        
    end
    
end

increase_summary = permute(increase_summary, [3 2 1]);

figure

for stat = 1:no_stats
    
    subplot(2, no_stats, stat)

    plot((1:no_blocks)', increase_summary(:, :, stat))
    
    hold on
    
    for block = 1:no_blocks
    
        plot(block, increase_summary(block, :, stat)', 'o')
        
    end
    
    set(gca, 'XLim', [.5 no_stats + .5], 'XTick', 1:no_stats, 'XTickLabel', {chan_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
    title(stat_labels{stat})
    
    subplot(2, no_stats, no_stats + stat)
    
    barwitherr(nanstd(increase_summary(:, :, stat)'), nanmean(increase_summary(:, :, stat)'));
    
    set(gca, 'XTickLabel', {chan_labels{:}, 'Overlap'}, 'FontSize', 14)
    
    ylabel('Time (min.)')
    
end

save([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix,...
    '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, norm,...
   win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_summary.mat'], 'increase_summary')

end
