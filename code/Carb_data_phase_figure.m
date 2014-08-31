function Carb_data_phase_figure

outlier_lim = 7; sd_lim = 2; win_size = 333; smooth_size = 20000;

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

f_bins = 9:2:31; no_f_bins = length(f_bins) - 1;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

f = 1000*(0:win_size)/win_size;

f_indices = f <= 32 & f >= 8;

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

record_label = {'st_m1', 'st_stn'};

figure;

%% Plotting phase angle by frequency, pre and post.

for r = 1:2
    
    load([record_label, '_subjects.mat'])
    
    load([record_label, '_', par_name, '_ch1_by_ch2_beta_ri_rose_dp.mat']);
    
    subplot(2, 2, r)
    
    h = barwitherr(conf_mat', angle(MR_mat)');
    
    bar_pos = get_bar_pos(h);
    
    f_bar_pairs = {};
    
    f_angle_indicator = nan(size(f_angle_pval));
    
    for pd = 1:2
        
        % Choose whichever is smaller - significant pairs or
        % insignificant pairs.
        if sum(f_angle_pval(:, pd) < 0.05) <= no_f_pairs/2
            
            f_angle_indicator(:, pd) = f_angle_pval(:, pd) < 0.05;
            
        else
            
            f_angle_indicator(:, pd) = f_angle_pval(:, pd) >= 0.05;
            
        end
        
        for fp = 1:no_f_pairs
            
            if f_angle_indicator(fp, pd)
                
                f_bar_pairs = {f_bar_pairs{:}, [bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 1)), bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 2))]};
                
            end
            
        end
        
    end
    
    f_angle_pval = reshape(f_angle_pval, 2*no_f_pairs, 1);
    
    f_angle_indicator = reshape(f_angle_indicator, 2*no_f_pairs, 1);
    
    sigstar(f_bar_pairs, f_angle_pval(f_angle_indicator == 1)')
    
    title(['Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'])
    
    colormap(c_order)
    
    colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
    
    set(gca, 'XTickLabel', period_label)
    
    %% Plotting coherence.
    
    subplot(2, 2, 2 + r)

    coh_listname = ['All_', record_label, '_ch', num2str(2 - r), '_', pd_label{pd}, '_coh_mtm_4tbw'];
    
    load(coh_listname)
    
    boundedline(f(f_indices)', mean_data(f_indices, :), std_data(f_indices, :, :))
    
    legend(period_label)
    
    axis tight
    
    xlabel('Frequency (Hz)')
    
    title({[chan_labels{2 - r}, ' High Beta Blocks'];'Phase of Coherence'})

end

save_as_pdf(gcf, 'Carb_phase_figure')

end

function pos_bars = get_bar_pos(handle)

for i = 1:length(handle)
    
    x = get(get(handle(i), 'children'), 'xdata');
    
    x = mean(x([1 3],:));
    
    pos_bars(i,:) = x;
    
end

end