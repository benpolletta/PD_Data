function Carb_data_phase_figure

outlier_lim = 7; sd_lim = 2; win_size = 333; smooth_size = 20000;

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

f_bins = 9:2:31; no_f_bins = length(f_bins) - 1;

f_pairs = nchoosek(1:no_f_bins, 2); no_f_pairs = size(f_pairs, 1);

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

f_center_indices = f_centers <= 25;

% c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

win_size = 2000;

f = 1000*(0:win_size)/win_size;

f_indices = f <= 26 & f >= 9;

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

record_label = {'st_m1', 'st_stn'}; record_chan_labels = {'_ch2_by_ch1_', '_ch1_by_ch2_'};

figure;

%% Plotting phase angle by frequency, pre and post.

for r = 1:2
    
    record_multiplier = (-1)^(r + 1);
    
    load([record_label{r}, '_subjects.mat'])
    
    load([record_label{r}, '_', par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat']);
    
    subplot(2, 2, r)
    
    conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
    
    conf_mat = repmat(conf_mat, [1 2 1]);
    
    h = boundedline(f_centers(f_center_indices)', record_multiplier*angle(MR_mat(f_center_indices, :)), conf_mat(f_center_indices, :, :));
    
    set(h, 'Marker', 's')
   
    axis tight
    
    hold on
    
    plot(f(f_indices)', zeros(length(f(f_indices)), 1), '--k')
    
    if r == 1
    
        legend(period_label, 'Location', 'SouthEast')
        
    end
    
    axis tight
    
    xlabel('Frequency (Hz)')
    
    % h = barwitherr(conf_mat', record_multiplier*angle(MR_mat)');
    % 
    % bar_pos = get_bar_pos(h);
    % 
    % f_bar_pairs = {};
    % 
    % f_angle_indicator = nan(size(f_angle_pval));
    % 
    % for pd = 1:2
    % 
    %     % Choose whichever is smaller - significant pairs or
    %     % insignificant pairs.
    %     if sum(f_angle_pval(:, pd) < 0.05) <= no_f_pairs/2
    % 
    %         f_angle_indicator(:, pd) = f_angle_pval(:, pd) < 0.05;
    % 
    %     else
    % 
    %         f_angle_indicator(:, pd) = f_angle_pval(:, pd) >= 0.05;
    % 
    %     end
    % 
    %     for fp = 1:no_f_pairs
    % 
    %         if f_angle_indicator(fp, pd)
    % 
    %             f_bar_pairs = {f_bar_pairs{:}, [bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 1)), bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 2))]};
    % 
    %         end
    % 
    %     end
    % 
    % end
    % 
    % f_angle_pval = reshape(f_angle_pval, 2*no_f_pairs, 1);
    % 
    % f_angle_indicator = reshape(f_angle_indicator, 2*no_f_pairs, 1);
    % 
    % sigstar(f_bar_pairs, f_angle_pval(f_angle_indicator == 1)')
    
    title({[chan_labels{3 - r}, ' High Beta Blocks']; ['Mean Phase Angle (', chan_labels{r}, ' - ', chan_labels{3 - r}, ') by ', chan_labels{r}, ' Freq.']})
    
    % colormap(c_order)
    % 
    % colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
    % 
    % set(gca, 'XTickLabel', period_label)
    
    %% Plotting coherence.
    
    subplot(2, 2, 2 + r)

    coh_listname = ['All_', record_label{r}, '_ch', num2str(3 - r), '_post_coh_mtm_4tbw_phase.mat'];
    
    load(coh_listname)
    
    boundedline(f(f_indices)', record_multiplier*mean_data(f_indices, :), std_data(f_indices, :, :))
    
    hold on
    
    plot(f(f_indices)', zeros(length(f(f_indices)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'NorthWest')
    
    end
        
    axis tight
    
    xlabel('Frequency (Hz)')
    
    title({[chan_labels{3 - r}, ' High Beta Blocks'];'Phase of Coherence'})

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