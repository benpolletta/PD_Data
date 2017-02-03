function power_post_vs_pre_timeseries_data(subject_mat_cell, file_name, peak_suffix, norm, freqs, no_cycles, bands, band_indices, window_length)

% Plots time series of beta power for all subjects.

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

if isempty(window_length) || window_length == 150
    
    window_length = 150; 

    win_flag = '';
    
else
    
    win_flag = ['_win', num2str(window_length)];
    
end

All_BP_plot = nan(30*60*sampling_freq, total_folders, no_chans);

[All_BP_increase, All_increase] = deal(nan(30*60*sampling_freq, total_folders, 2, no_chans)); 

BP_increased = zeros(total_folders, 2, no_chans);

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
            
            All_t_indices = All_time >= t(1) & All_time <= t(end);
            
            BP_t_indices = t >= All_time(1) & t <= All_time(end);
            
            load([subj_name, BP_suffix,win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_stats.mat'])
            
            blocks = cell(3, 1);
            
            for ch = 1:no_chans
                
                BP_plot = nanzscore(nanconv(BP(BP_t_indices, b, chan_order(ch)),...
                    ones(window_length*sampling_freq, 1)/(window_length*sampling_freq), 'nanout'));
                
                All_BP_plot(All_t_indices, folder_index, ch) = BP_plot;
                
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
                            
                            blocks{block_indices(bl)} = blocks{block_indices(bl)}';
                            
                        end
                        
                        increase_sec_blocks = blocks{block_indices(bl)} +...
                            (epoch_secs/2)*ones(size(blocks{block_indices(bl)}))*diag([-1 1]);
                        
                        increase_index = All_t_indices;
                        
                        increase_index(All_t_indices) = blocks_to_index(increase_sec_blocks, t(BP_t_indices));
                        
                        All_increase(:, folder_index, bl, ch) = increase_index;
                        
                        if sum(increase_index) == 0
                        
                            All_BP_increase(:, folder_index, bl, ch) = All_BP_plot(:, folder_index, chan_order(ch));
                        
                        else
                        
                            All_BP_increase(increase_index, folder_index, bl, ch) = All_BP_plot(increase_index, fo, chan_order(ch));
                        
                            BP_increased(folder_index, bl, ch) = 1;
                        
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        % subj_mat_name = subj_mat_cell{s};
        % 
        % subj_mat_name = subj_mat_name(1:(end - length('_subjects.mat')));
        
    end
    
    save([file_name, BP_suffix, win_flag, '_', short_band_labels{b} '_power_post_vs_pre_timeseries.mat'],...
        'All_time', 'All_BP_plot', 'All_increase', 'All_BP_increase', 'BP_increased')
    
end

end
