function PD_beta_epochs_rel_infusion_roseplot_by_dp_ind_avg(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

no_subs = length(folders);

f_bins = 8:4:32; no_f_bins = length(f_bins) - 1;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

load([subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_','_beta_ri_rose_dp.mat'])

mean_histogram = nanmean(All_histograms, 3);

mean_DP_mat = circ_mean(angle(All_MR_mat), [], 3);

for ch = 1:no_channels
    
    for ch1 = 1:2
        
        figure
        
        for pd = 1:2
            
            subplot(2, 2, pd)
            
            imagesc(bins{2}, bins{1}, mean_histogram(:, :, :, pd, ch, ch1))
            
            axis xy
            
            xlabel('Frequency (Hz)')
            
            ylabel('Phase Lag (rad)')
            
            title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Phase Lag by ', chan_labels{ch1}, ' Freq.']})
            
            freezeColors
            
        end
        
        mean_conf = nan(no_f_bins, 2);
        
        for f = 1:no_f_bins
            
            for pd = 1:2
                
                mean_conf(f, pd) = circ_confmean(reshape(angle(All_MR_mat(f, pd, :, ch, ch1)), size(All_MR_mat, 3), 1));
                
            end
            
        end
        
        mean_DP = reshape(mean_DP_mat(:, :, :, ch, ch1), no_f_bins, 2);
      
        c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';
        
        subplot(2, 3, 3 + 1)
        
        barwitherr(mean_conf, mean_DP)
        
        title(['Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'])
        
        set(gca, 'XTickLabel', f_centers)
        
        xlabel('Frequency')
        
        h = legend(period_label, 'Location', 'SouthEast'); freezeColors(h);
            
        freezeColors
        
        subplot(2, 3, 3 + 2)
        
        barwitherr(mean_conf', mean_DP')
       
        colormap(c_order)            
        
        title(['Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'])
        
        colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
        
        set(gca, 'XTick', [1 2], 'XTickLabel', period_label)
        
        subplot(2, 3, 3 + 3)
        
        plot(cumsum(ones(size(All_slopes, 1), 2), 2)', All_slopes(:, :, ch, ch1)', '*-')
        
        xlim([.8 2.2])
        
        title('Slope of Regression Line, Phase Angle vs. Freq.')
        
        set(gca, 'XTickLabel', period_label)
        
        save_as_pdf(gcf, [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_ri_rose_dp'])
        
    end
    
end

close('all')
    