function power_post_vs_pre_timeseries_plot(file_name, chan_labels, peak_suffix, freqs, no_cycles, bands, band_indices, window_length)

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

no_chans = length(chan_labels);

symbols = {'x', '+'};

included_folders = ones(14, 1);

included_folders([3 9]) = 0;

included_folders = logical(included_folders);

for b = band_indices
    
    name = [file_name, BP_suffix, win_flag, '_', short_band_labels{b} '_power_post_vs_pre_timeseries'];
    
    load([name, '.mat'])
    
    for ch = 1:no_chans
        
        figure(b)
        
        h = subplot(2, 1, ch);
        
        set(h, 'NextPlot', 'add', 'ColorOrder', distinguishable_colors(size(All_BP_plot, 2)), 'FontSize', 16)
        
        handle = plot(All_time/60, All_BP_plot(:, included_folders, ch), 'LineWidth', 1);
        
        axis tight
        
        box off
        
        change_line_style(handle, All_time/60, All_BP_plot(:, included_folders, ch),...
            All_increase(:, included_folders, 1, ch), 'linewidth', 2)
        
        % change_line_style(handle, All_time/60, All_BP_plot(:, included_folders, ch),...
        %    All_increase(:, included_folders, 2, ch), 'linewidth', 5)
        
        % hold on
        
        % add_stars_one_line(gca, All_time/60, All_increase(:, :, 1, ch), 1, 'symbol', symbols{ch})
        
        % for block = 1:2
        % 
        %     h = plot(All_time/60, All_BP_increase(:, :, block, ch));
        % 
        %     line_indices = find(BP_increased(:, block, ch));
        % 
        %     for l = 1:length(line_indices)
        % 
        %         set(h(line_indices(l)), 'LineWidth', 2*block)
        % 
        %     end
        % 
        % end
        
        if ch == 1
            
            title(sprintf('%d - %d Hz Power', bands(b, :)))
            
        elseif ch == no_chans
            
            xlabel('Time Rel. Infusion (Min.)')
            
        end
            
        ylabel({chan_labels{ch}; sprintf('Power per %g Min. (a.u.)', window_length/60)})
        
        plot([0; 0], [all_dimensions(@min, All_BP_plot); all_dimensions(@max, All_BP_plot)], 'k', 'LineWidth', 1)
        
        
    end
        
    save_as_pdf(gcf, name)
    
end

end
