function Carb_coh_figure(smooth_size, coh_win_size)

% Defaults.
if isempty(smooth_size), smooth_size = 5000; end
if isempty(coh_win_size), coh_win_size = 1000; end

outlier_lim = 7; sd_lim = 2; sampling_freq = 1000;

coh_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(coh_win_size),'win_',num2str(smooth_size),'smooth'];

f_coh = sampling_freq*(0:coh_win_size)/coh_win_size;

f_indices_coh = f_coh <= 30 & f_coh >= 10;

f_bins = [10 18 26]; no_bins = length(f_bins) - 1; f_bin_labels = cell(no_bins, 1);

for b = 1:no_bins
    
    f_bin_labels{b} = sprintf('%d - %d', f_bins(b), f_bins(b + 1));
    
end

period_label = {'Pre-Infusion', 'Post-Infusion'};

record_label = {'st_m1', 'st_stn'};
    
rad_deg = 180/pi;

fig_par_name = [num2str(smooth_size),'smooth_',num2str(coh_win_size),'coh'];

figure;

for r = 1:2
    
    record_multiplier = (-1)^(r + 1);
    
    load([record_label{r}, '_subjects.mat'])
    
    %% Plotting coherence.
    
    subplot(2, 2, r)

    coh_listname = ['All_', record_label{r}, '_ch', num2str(3 - r), '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'];
    
    load(coh_listname)
    
    h = boundedline(f_coh(f_indices_coh)', rad_deg*record_multiplier*mean_data(f_indices_coh, :), rad_deg*std_data(f_indices_coh, :, :));
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'SouthEast')
    
    end
        
    axis tight
    
    title({[chan_labels{3 - r}, ' High Beta Blocks'];'Phase of Coherence'})

    subplot(2, 2, 2 + r)
    
    [bin_mean_coh, bin_std_coh] = deal(nan(no_bins, size(mean_data, 2)));
    
    for b = 1:no_bins
        
        bin_indices = f_bins(b) <= f_coh & f_coh <= f_bins(b + 1);
        
        bin_mean_coh(b, :) = record_multiplier*nanmean(mean_data(bin_indices, :));
        
        std_sum = nansum(std_data(bin_indices, :, 1).^2);
        
        mean_sum = nansum((mean_data(bin_indices, :) - ones(size(mean_data(bin_indices, :)))*diag(bin_mean_coh(b, :))).^2);
        
        bin_std_coh(b, :) = sqrt((std_sum + mean_sum)/length(bin_indices));
        
    end
    
    h = barwitherr(bin_std_coh, bin_mean_coh);
    
    set(h(2), 'FaceColor', [0 .5 0])
    
    set(gca, 'XTickLabel', f_bin_labels)
    
    xlabel('Frequency (Hz)')
    
    ylabel({['GC, (', chan_labels{r}, '-->', chan_labels{3 - r}, ')'];['- (', chan_labels{3 - r}, '-->', chan_labels{r}, ')']})

end

save_as_pdf(gcf, ['Carb_coh_figure_', fig_par_name])

end