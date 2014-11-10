function Carb_rp_figure(smooth_size, rp_win_size)

% Defaults.
if isempty(smooth_size), smooth_size = 5000; end
if isempty(rp_win_size), rp_win_size = 333; end

outlier_lim = 7; sd_lim = 2;

rp_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(rp_win_size),'win_',num2str(smooth_size),'smooth'];

f_bins{1} = 9.5:1:30.5; f_bins{2} = [10 18 26];

[f_centers, f_center_indices, no_bins] = deal(cell(2, 1));

for b = 1:2
    
    f_centers{b} = (f_bins{b}(1:(end-1)) + f_bins{b}(2:end))/2;
    
    f_center_indices{b} = f_centers{b} <= 30;
    
    no_bins{b} = length(f_center_indices{b});
    
end

period_label = {'Pre-Infusion', 'Post-Infusion'};

record_label = {'st_m1', 'st_stn'}; record_chan_labels = {'_ch2_by_ch1_', '_ch1_by_ch2_'};
    
rad_deg = 180/pi;

fig_par_name = [num2str(smooth_size),'smooth_',num2str(rp_win_size),'rp'];

figure;

for r = 1:2
    
    for b = 1:2
        
        record_multiplier = (-1)^(r + 1);
        
        load([record_label{r}, '_subjects.mat'])
        
        %% Plotting phase angle by frequency, pre and post.
        
        load([record_label{r}, '_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_', num2str(no_bins{b}), 'bins_group.mat'])
        
        subplot(2, 2, (b - 1)*2 + r)
        
        if b == 1
        
        conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
        
        conf_mat = repmat(conf_mat, [1 2 1]);
        
        h = boundedline(f_centers{b}(f_center_indices{b})', rad_deg*record_multiplier*angle(MR_mat(f_center_indices{b}, :)), rad_deg*conf_mat(f_center_indices{b}, :, :));
        
        set(h, 'Marker', 's')
        
        axis tight
        
        hold on
        
        plot(f_centers{b}(f_center_indices{b})', zeros(length(f_centers{b}(f_center_indices{b})), 1), '--k')
        
        else
           
            h = barwitherr(rad_deg*conf_mat, rad_deg*record_multiplier*angle(MR_mat));
            
            set(h(2), 'FaceColor', [0 .5 0])
            
            set(gca, 'XTickLabel', {'10 - 18', '18 - 26'})
            
        end
        
        if r == 1
            
            legend(period_label, 'Location', 'SouthEast')
            
        end
        
        title([chan_labels{3 - r}, ' High Beta Blocks'])
        
        ylabel(['Mean Phase Angle (', chan_labels{r}, ' - ', chan_labels{3 - r}, ')'])
        
        xlabel([chan_labels{r}, ' Freq.'])
        
    end
    
end

save_as_pdf(gcf, ['Carb_rp_figure_', fig_par_name])

end