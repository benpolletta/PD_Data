function Carb_data_lag_figure(smooth_size, rp_win_size, coh_win_size)

if isempty(smooth_size), smooth_size = 20000; end
if isempty(rp_win_size), rp_win_size = 333; end
if isempty(coh_win_size), coh_win_size = 2000; end

outlier_lim = 7; sd_lim = 2;

rp_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(rp_win_size),'win_',num2str(smooth_size),'smooth'];

f_bins = 9.5:1:30.5;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

rp_phase_time = (1000./f_centers)/(2*pi);

f_center_indices = f_centers <= 30;

coh_par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(coh_win_size),'win_',num2str(smooth_size),'smooth'];

f_coh = 1000*(0:coh_win_size)/coh_win_size;

coh_phase_time = (1000./f_coh)/(2*pi);

f_indices_coh = f_coh <= 30 & f_coh >= 10;

period_label = {'Pre-Infusion','Post-Infusion','6OHDA'};

record_label = {'st_m1', 'st_stn', 'st_m1_6OHDA'}; record_chan_labels = {'_ch2_by_ch1_', '_ch1_by_ch2_', '_ch2_by_ch1_'};

record_indices = {1:2, 1:2, 1}; record_channels = {'2', '1', '2'};

fig_legend = {'Str. - Motor Lag, Motor High \beta, Pre-Infusion', 'Str. - Motor Lag, Motor High \beta, Post-Infusion', 'Str. - STN Lag, STN High \beta, Pre-Infusion', 'Str. - STN Lag, STN High \beta, Post-Infusion', 'Str. - Motor Lag, Motor High \beta, 6OHDA'};

color_map = [0 0 1; 0 .5 0; .75 0 .75; 0 .75 .75; 1 0 0];

fig_par_name = [num2str(smooth_size), 'smooth_', num2str(rp_win_size), 'rp_', num2str(coh_win_size), 'coh'];

figure;

[All_MR_mat, All_conf_mat, All_mean_data, All_std_data] = deal([]);

for r = 1:3
    
    record_multiplier = (-1)^(r + 1);

    %% Collecting phase angle by frequency, pre and post.
    
    load([record_label{r}, '_subjects.mat'])
    
    load([record_label{r}, '_', rp_par_name, record_chan_labels{r}, 'beta_ri_rose_dp_group.mat'])
        
    All_MR_mat(:, (end + 1):(end + length(record_indices{r}))) = record_multiplier*angle(MR_mat(:, record_indices{r})); 
    
    All_conf_mat(:, (end + 1):(end + length(record_indices{r}))) = conf_mat(:, record_indices{r});
    
    %% Collecting coherence.

    coh_listname = ['All_', record_label{r}, '_ch', record_channels{r}, '_post_', coh_par_name, '_coh_mtm_', num2str(2*coh_win_size/1000), 'tbw_phase.mat'];
    
    load(coh_listname)
    
    All_mean_data(:, (end + 1):(end + length(record_indices{r}))) = record_multiplier*mean_data(:, record_indices{r});
    
    All_std_data(:, :, (end + 1):(end + length(record_indices{r}))) = std_data(:, :, record_indices{r});
    
end

All_std_data(:, :, 1) = [];

lag_const = [0, 2*pi];

lag_title = {'Motor Ctx. Leads at Short Lag', 'Striatum Leads at Long Lag'};

for l = 1:2
    
    %% Plotting lags from roseplots.
    
    subplot(2, 2, l)
    
    conf_plot = diag(rp_phase_time(f_center_indices))*All_conf_mat(f_center_indices, :);
    
    conf_plot = reshape(conf_plot, size(conf_plot, 1), 1, size(conf_plot, 2));
    
    conf_plot = repmat(conf_plot, [1 2 1]);
    
    h = boundedline(f_centers(f_center_indices)', diag(rp_phase_time(f_center_indices))*(lag_const(l) + All_MR_mat(f_center_indices, :)), conf_plot, 'cmap', color_map);
    
    set(h, 'Marker', 's')
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
    
    if l == 1
    
        legend(fig_legend, 'Location', 'SouthEast')
        
    end
    
    axis tight
    
    title({'Time Lag From Rose Plot,'; lag_title{l}; [num2str(smooth_size/1000), ' Smoothing']})
    
    %% Plotting lags from coherence.
    
    subplot(2, 2, 2 + l)
    
    std_plot = nan(sum(f_indices_coh), size(All_std_data, 2));
    
    for d3 = 1:size(All_std_data, 3)
       
        std_plot(:, :, d3) = diag(coh_phase_time(f_indices_coh))*ones(size(All_std_data(f_indices_coh, :, d3))).*All_std_data(f_indices_coh, :, d3);
        
    end
    
    h = boundedline(f_coh(f_indices_coh)', diag(coh_phase_time(f_indices_coh))*(lag_const(l) + All_mean_data(f_indices_coh, :)), std_plot, 'cmap', color_map);
    
    set(h, 'Marker', 's')
    
    hold on
    
    plot(f_coh(f_indices_coh)', zeros(length(f_coh(f_indices_coh)), 1), '--k')
        
    axis tight
    
    xlabel('Frequency (Hz)')
    
    title({'Time Lag From Coherence,'; lag_title{l}})
    
end

save_as_pdf(gcf, ['Carb_lag_figure_', fig_par_name])

end