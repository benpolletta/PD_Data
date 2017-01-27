function power_post_vs_pre_plot(subject_mat_cell, peak_suffix, epoch_secs, pd_handle, norm, freqs, no_cycles, bands, band_indices, window_length)

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200; in_freqs = [];
    
    no_cycles = linspace(3, 21, length(freqs)); in_no_cycles = [];
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; in_bands = [];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    in_freqs = freqs; in_no_cycles = no_cycles; in_bands = bands;
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

for s = 1:length(subject_mat_cell)
    
    load(subject_mat_cell{s})
    
    no_folders(s) = length(folders);
    
end

total_folders = sum(no_folders);

sampling_freq = 500; % load([folders{1}, '/', prefixes{1}, BP_suffix, '_wt.mat'], 'sampling_freq')

no_bands = size(bands, 1);

[short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

no_chans = length(chan_labels);

no_pds = length(pd_labels);

if isempty(window_length)
    
    window_length = epoch_secs; 

    win_flag = '';
    
else
    
    win_flag = ['_win', num2str(window_length)];
    
end

All_BP_plot = nan(30*60*sampling_freq, total_folders, 2);

All_BP_increase = nan(30*60*sampling_freq, total_folders, 2, 2); 

All_time = (1:(30*60*sampling_freq))/sampling_freq - 10*60;

for b = band_indices
    
    folder_index = 0;
    
    for s = 1:length(subject_mat_cell)
        
        load(subject_mat_cell{s})
        
        chan_order = [1 2];
        
        if ~strcmp(chan_labels{1}, 'Striatum');
            
            chan_order = fliplr(chan_order);
            
        end
        
        for fo = 1:no_folders(s)
            
            folder_index = folder_index + 1;
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            outlier_lim = outlier_lims(fo);
            
            subj_name = [folder,'/',prefix];
            
            if strcmp(norm, '_peaks')
                
                data = load([subj_name, BP_suffix, '_wt_BP.mat']);
                
                Spike_indicator = peak_loader(subj_name, peak_suffix, size(data.BP, 1));
                
                BP = repmat(permute(Spike_indicator, [1 2 3]), [1 no_bands 1]);
                
            else
                
                BP = get_BP(subj_name, peak_suffix, [], outlier_lim, norm, in_freqs, in_no_cycles, in_bands);
                
                % load([subj_name, BP_suffix, '_wt_BP.mat'])
                %
                % eval(['BP = BP', norm, ';'])
                
            end
            
            t = (1:size(BP, 1))/sampling_freq - basetimes(fo); % t = t/60;
            
            t_indices = All_time >= t(1) & All_time <= t(end);
            
            load([subj_name, BP_suffix,win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_stats.mat'])
            
            blocks = cell(3, 1);
            
            for ch = 1:no_chans
                
                BP_plot = nanzscore(nanconv(BP(:, b, chan_order(ch)), ones(5*60*sampling_freq, 1)/(5*60*sampling_freq), 'nanout'));
                
                All_BP_plot(t_indices, folder_index, ch) = BP_plot;
                
                increase_blocks = index_to_blocks(test_win(:, 2, chan_order(ch)));
                
                channel_sec_blocks = win_secs(increase_blocks);
                
                blocks{ch} = channel_sec_blocks;
                
            end
            
            overlap_index = test_win(:, 2, 1) & test_win(:, 2, 2);
            
            overlap_blocks = index_to_blocks(overlap_index);
            
            overlap_sec_blocks = win_secs(overlap_blocks);
            
            blocks{3} = overlap_sec_blocks;
            
            for ch = 1:no_chans
                
                block_indices = [ch 3];
                
                for bl = 1:2;
                    
                    if ~isempty(blocks{block_indices(bl)})
                        
                        if size(blocks{block_indices(bl)}, 2) ~= 2
                            
                            blocks{bl} = blocks{block_indices(bl)}';
                            
                        end
                        
                        increase_sec_blocks = blocks{bl} + (epoch_secs/2)*ones(size(blocks{bl}))*diag([-1 1]);
                        
                        increase_index = t_indices;
                        
                        increase_index(t_indices) = blocks_to_index(increase_sec_blocks, t);
                        
                        All_BP_increase(increase_index, folder_index, bl, ch) = All_BP_plot(increase_index, fo, ch);
                        
                    end
                    
                end
                
            end
            
        end
        
        for ch = 1:no_chans
            
            figure(b)
            
            subplot(2, 1, ch);
            
            plot(All_time/60, All_BP_plot(:, :, ch))
            
            axis tight
            
            box off
            
            hold on
            
            plot(All_time/60, All_BP_increase(:, :, 1, ch), 'LineWidth', 2)
            
            plot(All_time/60, All_BP_increase(:, :, 2, ch), 'LineWidth', 3)
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
                
            elseif fo == no_folders
                
                xlabel('Time Rel. Infusion (Min.)')
                
            end
            
            if ch == 1
                
                ylabel({folder; 'Power (per 5 Min.)'})
                
            end
            
            plot([0; 0], [min(BP_plot); max(BP_plot)], 'k', 'LineWidth', 1)
            
            
        end
        
    end
        
    save_as_pdf(gcf, [subj_name, BP_suffix, win_flag, '_', short_band_labels{b} '_power_post_vs_pre_timeseries.mat'])
    
end

end
