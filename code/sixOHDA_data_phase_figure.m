function sixOHDA_data_phase_figure(smooth_size, rp_win_size, coh_win_size, GC_win_size)

% Defaults.
if isempty(smooth_size), smooth_size = 5000; end
if isempty(rp_win_size), rp_win_size = 333; end
if isempty(coh_win_size), coh_win_size = 1000; end
if isempty(GC_win_size), GC_win_size = 1000; end

outlier_lim = 7; sd_lim = 2; sampling_freq = 1000;

rp_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(rp_win_size),'win_',num2str(smooth_size),'smooth'];

f_bins = 9.5:1:30.5;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

f_center_indices = f_centers <= 30;

coh_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(coh_win_size),'win_',num2str(smooth_size),'smooth'];

f_coh = sampling_freq*(0:coh_win_size)/coh_win_size;

f_indices_coh = f_coh <= 30 & f_coh >= 10;

GC_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(GC_win_size),'win_',num2str(smooth_size),'smooth'];

f_GC = sampling_freq*(0:GC_win_size)/GC_win_size;

f_indices_GC = f_GC <= 30 & f_GC >= 10;

period_label = {'Pre-Infusion', '6OHDA'};

record_label = {'st_m1', 'st_stn'}; record_chan_labels = {'_ch2_by_ch1_', '_ch1_by_ch2_'};
    
rad_deg = 180/pi;

fig_par_name = [num2str(smooth_size),'smooth_',num2str(rp_win_size),'rp_',num2str(coh_win_size),'coh_',num2str(GC_win_size),'GC'];

color_map = [0 0 1; 1 0 0];

figure;

for r = 1:1
    
    record_multiplier = (-1)^(r + 1);

    %% Plotting phase angle by frequency, pre and post.
    
    load([record_label{r}, '_subjects.mat'])
    
    load([record_label{r}, '_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat'])
    
    if r == 1
        
        MR_mat_temp = MR_mat; conf_mat_temp = conf_mat;
        
        load([record_label{r}, '_6OHDA_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat'])
        
        MR_mat_temp(:, 2) = MR_mat(:, 1); conf_mat_temp(:, 2) = conf_mat(:, 1);
        
        clear MR_mat conf_mat
        
        MR_mat = MR_mat_temp; conf_mat = conf_mat_temp;
        
    end
    
    subplot(3, 1, r)
    
    conf_mat = reshape(conf_mat, size(conf_mat, 1), 1, size(conf_mat, 2));
    
    conf_mat = repmat(conf_mat, [1 2 1]);
    
    h = boundedline(f_centers(f_center_indices)', rad_deg*record_multiplier*angle(MR_mat(f_center_indices, :)), rad_deg*conf_mat(f_center_indices, :, :), 'cmap', color_map);
    
    set(h, 'Marker', 's')
   
    axis tight
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if r == 1
    
        legend(period_label, 'Location', 'SouthEast')
        
    end
    
    axis tight
    
    title({[chan_labels{3 - r}, ' High Beta Blocks']; ['Mean Phase Angle (', chan_labels{r}, ' - ', chan_labels{3 - r}, ') by ', chan_labels{r}, ' Freq.']})
    
    %% Plotting coherence.
    
    subplot(3, 1, 1 + r)

    coh_listname = ['All_', record_label{r}, '_ch', num2str(3 - r), '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'];
    
    load(coh_listname)
    
    if r == 1
        
        mean_data_temp = mean_data; std_data_temp = std_data;
        
        load(['All_', record_label{r}, '_6OHDA', '_ch', num2str(3 - r), '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'])
        
        mean_data_temp(:, 2) = mean_data(:, 1); std_data_temp(:, :, 2) = std_data(:, :, 1);
        
        clear mean_data std_data
        
        mean_data = mean_data_temp; std_data = std_data_temp;
        
    end
    
    h = boundedline(f_coh(f_indices_coh)', rad_deg*record_multiplier*mean_data(f_indices_coh, :), rad_deg*std_data(f_indices_coh, :, :), 'cmap', color_map);
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if r == 1
   
        legend(period_label, 'Location', 'SouthEast')
    
    end
        
    axis tight
    
    title({[chan_labels{3 - r}, ' High Beta Blocks'];'Phase of Coherence'})
    
    %% Plotting Granger causality.
    
    subplot(3, 1, 2 + r)
    
    load(['All_', record_label{r}, '_', GC_par_name, '_GC_spec_stats.mat'])
    
    plot_mean_GC = All_mean_GC_spec(:, :, 2 - r, 1);
    
    plot_std_GC = All_std_GC_spec(:, :, :, 2 - r, 1);
    
    if r == 1
    
        plot_mean_GC_temp = plot_mean_GC; plot_std_GC_temp = plot_std_GC;
    
        load(['All_', record_label{r}, '_6OHDA_', GC_par_name, '_GC_spec_stats.mat'])
    
        plot_mean_GC_temp(:, 2) = All_mean_GC_spec(:, 1, 3 - r, 1); plot_std_GC_temp(:, :, 2) = All_std_GC_spec(:, :, 1, 3 - r, 1);
    
        clear plot_mean_GC plot_std_GC
    
        plot_mean_GC = plot_mean_GC_temp; plot_std_GC = plot_std_GC_temp;
    
    end
    
    h = boundedline(f_GC(f_indices_GC)', record_multiplier*plot_mean_GC(f_indices_GC, :), plot_std_GC(f_indices_GC, :, :), 'cmap', color_map);
    
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

save_as_pdf(gcf, ['sixOHDA_phase_figure_', fig_par_name])

end

% function pos_bars = get_bar_pos(handle)
% 
% for i = 1:length(handle)
%     
%     x = get(get(handle(i), 'children'), 'xdata');
%     
%     x = mean(x([1 3],:));
%     
%     pos_bars(i,:) = x;
%     
% end
% 
% end