function power_post_vs_pre_stats(subject_mat, peak_suffix, epoch_secs, pd_handle, norm, in_freqs, in_no_cycles, in_bands, band_indices, window_length)

[~, ~, bands, BP_suffix] = init_freqs(in_freqs, in_no_cycles, in_bands);

BP_suffix = [BP_suffix, peak_suffix];

load(subject_mat)

no_folders = length(folders);

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

window_step = min(30, window_length);

load([subject_mat(1:(end - length('_subjects.mat'))), BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs', pd_handle, '.mat'])

for fo = 1:no_folders
    
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
    
    clear pd_indices
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0; pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), 1);
        
    end
    
    pd_length = nan(no_pds, 1);
    
    for pd = 1:no_pds
        
        pd_length(pd) = sum(pd_indices(:, pd));
        
        if pd == 1
            
            display(sprintf('Baseline recording length for %s = %g minutes.', folder, pd_length(1)/(60*sampling_freq)))

        end
        
    end
    
    if no_pds == 2
    
        no_post_windows = floor((pd_length(2) - (window_length*sampling_freq)/2)/(window_step*sampling_freq));
        
    else
        
        no_post_windows = 0;
        
    end
    
    pd_indices = logical(pd_indices);
    
    for b = band_indices
            
        [p_win, test_win] = deal(nan(no_post_windows, 2, no_chans));
        
        for ch = 1:no_chans
            
            handle = nan(no_pds, 1);
            
            BP_plot = nanconv(BP(:, b, ch), ones(5*60*sampling_freq, 1)/(5*60*sampling_freq), 'nanout');
            
            [BP_sec, t_sec] = deal(cell(2, 1));
            
            for pd = 1:no_pds
                
                pd_secs = floor(pd_length(pd)/sampling_freq);
                
                [BP_sec{pd}, t_sec{pd}] = deal(nan(pd_secs, 1));
                
                BP_selected = BP(pd_indices(:, pd), b, ch);
                
                if pd == 1
                    
                    BP_selected = flipud(BP_selected);
                    
                end
                
                BP_selected = nans_to_end(BP_selected);
                
                for sec = 1:pd_secs
                    
                    sec_start = (sec - 1)*sampling_freq + 1;
                    
                    sec_end = min(sec*sampling_freq, size(BP_selected, 1));
                    
                    t_sec{pd}(sec) = (sec_start + sec_end)/(2*sampling_freq);
                    
                    BP_sec{pd}(sec) = nanmean(BP_selected(sec_start:sec_end));
                    
                end
                
                if pd == 1
                    
                    t_sec{pd} = -t_sec{pd};
                    
                end
                
            end
            
            [win_centers, win_secs] = deal(nan(no_post_windows, 1));
            
            for w = 1:no_post_windows
                
                win_start = (w - 1)*window_step + 1;
                
                win_end = (w - 1)*window_step + window_length;
                
                win_centers(w) = (win_start + win_end)/2;
                
                win_indices = t_sec{2} >= win_start & t_sec{2} <= win_end;
                
                win_secs(w) = nanmean(t_sec{2}(win_indices));
                
                [~, p_win(w, 1, ch)] = ttest2(BP_sec{1}, BP_sec{2}(win_indices), 'tail', 'right');
                
                [~, p_win(w, 2, ch)] = ttest2(BP_sec{1}, BP_sec{2}(win_indices), 'tail', 'left');
                
            end
            
            test_win(:, :, ch) = p_win(:, :, ch) < .05;
        
            figure(b)
            
            h = subplot(no_folders, 2, (fo - 1)*2 + ch);
            
            plot(t/60, BP_plot)
            
            axis tight
            
            box off
            
            hold on
            
            add_stars(h, win_centers/60, test_win(:, :, ch), [0 1], [1 0 0; 1 .5 0])
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
                
            elseif fo == no_folders
                
                xlabel('Time Rel. Infusion (Min.)')
                
            end
            
            if ch == 1
                
                ylabel({folder; 'Power (per Min.)'})
                
            end
            
            plot([0; 0], [min(BP_plot); max(BP_plot)], 'k', 'LineWidth', 1)
            
            figure(no_bands + b)
            
            subplot(no_folders, 2, (fo - 1)*2 + ch)
            
            semilogy(win_centers/60, p_win(:, :, ch))
            
            axis tight
            
            box off
            
            hold on
            
            semilogy(range(win_centers/60)', .05*ones(2, 1), 'k')
            
            if fo == 1
                
                title(sprintf('%s, %d - %d Hz', chan_labels{ch}, bands(b, :)))
                
            elseif fo == no_folders
                
                xlabel('Window Center (Min. Rel. Infusion)')
                
            end
            
            if ch == 1
                
                ylabel({folder; 'p-Value, Window vs. Baseline'})
                
            end
            
            plot([0; 0], [min(BP_plot); max(BP_plot)], 'k', 'LineWidth', 1)
            
        end
            
        save([subj_name, BP_suffix, win_flag, '_', short_band_labels{b}, '_power_post_vs_pre_stats.mat'], 't_sec', 'p_win', 'win_secs', 'test_win')
        
    end
    
end

for b = band_indices
           
    save_as_pdf(b, [subject_mat(1:(end - length('_subjects.mat'))), '_', short_band_labels{b},...
        pd_handle, norm, win_flag, '_power_post_vs_pre'])
    
    save_as_pdf(no_bands + b, [subject_mat(1:(end - length('_subjects.mat'))), '_', short_band_labels{b},...
        pd_handle, norm, win_flag, '_power_post_vs_pre_pvals'])
    
end

end
