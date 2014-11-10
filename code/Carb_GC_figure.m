function Carb_GC_figure(smooth_size, GC_win_size)

% Defaults.
if isempty(smooth_size), smooth_size = 5000; end
if isempty(GC_win_size), GC_win_size = 1000; end

outlier_lim = 7; sd_lim = 2; sampling_freq = 1000;

GC_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(GC_win_size),'win_',num2str(smooth_size),'smooth'];

f_GC = sampling_freq*(0:GC_win_size)/GC_win_size;

f_indices_GC = f_GC <= 31 & f_GC >= 9;

f_bins = [10 18 26]; no_bins = length(f_bins) - 1; f_bin_labels = cell(no_bins, 1);

for b = 1:no_bins
    
    f_bin_labels{b} = sprintf('%d - %d', f_bins(b), f_bins(b + 1));
    
end

period_label = {'Pre-Infusion', 'Post-Infusion'};

record_label = {'st_m1', 'st_stn'};

fig_par_name = [num2str(smooth_size),'smooth_',num2str(GC_win_size),'GC'];

figure;

for r = 1:2
    
    record_multiplier = (-1)^(r + 1);
    
    load([record_label{r}, '_subjects.mat'])
    
    %% Plotting Granger causality.
    
    subplot(2, 2, r)%(3, 2, 4 + r)
    
    load(['All_', record_label{r}, '_', GC_par_name, '_GC_spec_stats.mat'])
    
    plot_mean_GC = All_mean_GC_spec(:, :, 3 - r, 1);
    
    plot_std_GC = All_std_GC_spec(:, :, :, 3 - r, 1);
    
    h = boundedline(f_GC(f_indices_GC)', record_multiplier*plot_mean_GC(f_indices_GC, :), plot_std_GC(f_indices_GC, :, :));
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_GC(f_indices_GC)', zeros(length(f_GC(f_indices_GC)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'SouthEast')
    
    end
        
    axis tight
    
    %xlabel('Frequency (Hz)')
    
    title([chan_labels{3 - r}, ' High Beta Blocks'])
    
    ylabel({['GC, (', chan_labels{r}, '-->', chan_labels{3 - r}, ')'],['- (', chan_labels{3 - r}, '-->', chan_labels{r}, ')']})

    subplot(2, 2, 2 + r)
    
    [bin_mean_GC, bin_std_GC] = deal(nan(no_bins, size(plot_mean_GC, 2)));
    
    for b = 1:no_bins
        
        bin_indices = f_bins(b) <= f_GC & f_GC <= f_bins(b + 1);
        
        bin_mean_GC(b, :) = record_multiplier*nanmean(plot_mean_GC(bin_indices, :));
        
        std_sum = nansum(plot_std_GC(bin_indices, :, 1).^2);
        
        mean_sum = nansum((plot_mean_GC(bin_indices, :) - ones(size(plot_mean_GC(bin_indices, :)))*diag(bin_mean_GC(b, :))).^2);
        
        bin_std_GC(b, :) = sqrt((std_sum + mean_sum)/length(bin_indices));
        
    end
    
    h = barwitherr(bin_std_GC, bin_mean_GC);
    
    set(h(2), 'FaceColor', [0 .5 0])
    
    set(gca, 'XTickLabel', f_bin_labels)
    
    xlabel('Frequency (Hz)')
    
    ylabel({['GC, (', chan_labels{r}, '-->', chan_labels{3 - r}, ')'];['- (', chan_labels{3 - r}, '-->', chan_labels{r}, ')']})
    
end

save_as_pdf(gcf, ['Carb_GC_figure_', fig_par_name])

end