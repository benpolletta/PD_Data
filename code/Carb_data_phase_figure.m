function Carb_data_phase_figure

outlier_lim = 7; sd_lim = 2; rp_win_size = 333; smooth_size = 20000;

rp_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(rp_win_size),'win_',num2str(smooth_size),'smooth'];

sixOHDA_rp_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(rp_win_size),'win_',num2str(5000),'smooth'];

f_bins = 9.5:1:30.5; no_f_bins = length(f_bins) - 1;

f_pairs = nchoosek(1:no_f_bins, 2); no_f_pairs = size(f_pairs, 1);

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

f_center_indices = f_centers <= 30;

% c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

coh_win_size = 1000;

coh_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(coh_win_size),'win_',num2str(smooth_size),'smooth'];

sixOHDA_coh_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(coh_win_size),'win_',num2str(5000),'smooth'];

f_coh = 1000*(0:coh_win_size)/coh_win_size;

f_indices_coh = f_coh <= 30 & f_coh >= 10;

f_GC = 100*(0:rp_win_size)/rp_win_size;

f_indices_GC = f_GC <= 30 & f_GC >= 10;

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion','6OHDA'};

record_label = {'st_m1', 'st_stn'}; record_chan_labels = {'_ch2_by_ch1_', '_ch1_by_ch2_'};
    
rad_deg = 180/pi;

figure;

for r = 1:2
    
    record_multiplier = (-1)^(r + 1);

    %% Plotting phase angle by frequency, pre and post.
    
    load([record_label{r}, '_subjects.mat'])
    
    load([record_label{r}, '_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat'])
    
    if r == 1
        
        MR_mat_temp = MR_mat; conf_mat_temp = conf_mat;
        
        load([record_label{r}, '_6OHDA_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat'])
        
        MR_mat_temp(:, 3) = MR_mat(:, 1); conf_mat_temp(:, 3) = conf_mat(:, 1);
        
        clear MR_mat conf_mat
        
        MR_mat = MR_mat_temp; conf_mat = conf_mat_temp;
        
    end
    
    subplot(3, 2, r)
    
    conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
    
    conf_mat = repmat(conf_mat, [1 2 1]);
    
    h = boundedline(f_centers(f_center_indices)', rad_deg*record_multiplier*angle(MR_mat(f_center_indices, :)), rad_deg*conf_mat(f_center_indices, :, :));
    
    set(h, 'Marker', 's')
   
    axis tight
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if r == 1
    
        legend(period_label, 'Location', 'SouthEast')
        
    end
    
    axis tight
    
    %xlabel('Frequency (Hz)')
    
    title({[chan_labels{3 - r}, ' High Beta Blocks']; ['Mean Phase Angle (', chan_labels{r}, ' - ', chan_labels{3 - r}, ') by ', chan_labels{r}, ' Freq.']})
    
    % colormap(c_order)
    % 
    % colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
    % 
    % set(gca, 'XTickLabel', period_label)
    
    %% Plotting coherence.
    
    subplot(3, 2, 2 + r)

    coh_listname = ['All_', record_label{r}, '_ch', num2str(3 - r), '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'];
    
    load(coh_listname)
    
    if r == 1
        
        mean_data_temp = mean_data; std_data_temp = std_data;
        
        load(['All_', record_label{r}, '_6OHDA', '_ch', num2str(3 - r), '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'])
        
        mean_data_temp(:, 3) = mean_data(:, 1); std_data_temp(:, :, 3) = std_data(:, :, 1);
        
        clear mean_data std_data
        
        mean_data = mean_data_temp; std_data = std_data_temp;
        
    end
    
    h = boundedline(f_coh(f_indices_coh)', rad_deg*record_multiplier*mean_data(f_indices_coh, :), rad_deg*std_data(f_indices_coh, :, :));
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'SouthEast')
    
    end
        
    axis tight
    
    %xlabel('Frequency (Hz)')
    
    title({[chan_labels{3 - r}, ' High Beta Blocks'];'Phase of Coherence'})
    
    %% Plotting Granger causality.
    
    subplot(3, 2, 4 + r)
    
    load(['All_', record_label{r}, '_', rp_par_name, '_GC_spec_stats.mat'])
    
    plot_mean_GC = All_mean_GC_spec(:, :, 3 - r, 1);
    
    plot_std_GC = All_std_GC_spec(:, :, :, 3 - r, 1);
    
    % if r == 1
    % 
    %     plot_mean_GC_temp = plot_mean_GC; plot_std_GC_temp = plot_std_GC;
    % 
    %     load(['All_', record_label{r}, '_60HDA_', rp_par_name, '_GC_spec_stats.mat'])
    % 
    %     plot_mean_GC_temp(:, 3) = All_mean_GC_spec(:, 1, 3 - r, 1); std_GC_temp(:, :, 3) = All_std_GC_spec(:, :, 1, 3 - r, 1);
    % 
    %     clear plot_mean_GC plot_std_GC
    % 
    %     plot_mean_GC = plot_mean_GC_temp; plot_std_GC = plot_std_GC_temp;
    % 
    % end
    
    h = boundedline(f_GC(f_indices_GC)', record_multiplier*plot_mean_GC(f_indices_GC, :), plot_std_GC(f_indices_GC, :, :));
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_GC(f_indices_GC)', zeros(length(f_GC(f_indices_GC)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'SouthEast')
    
    end
        
    axis tight
    
    xlabel('Frequency (Hz)')
    
    title({[chan_labels{3 - r}, ' High Beta Blocks'];['GC, (', chan_labels{r}, '-->', chan_labels{3 - r}, ') - (', chan_labels{3 - r}, '-->', chan_labels{r}, ')']})

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