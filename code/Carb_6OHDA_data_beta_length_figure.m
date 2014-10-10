function Carb_6OHDA_data_beta_length_figure(win_size, smooth_size, no_bins)

if isempty(smooth_size), smooth_size = 20000; end
if isempty(win_size), win_size = 333; end

outlier_lim = 7; sd_lim = 2;

if isscalar(win_size)

    par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

    fig_par_name = [num2str(smooth_size), 'smooth_', num2str(win_size), 'win'];

elseif length(win_size) == 2
   
    par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size(1)), 'to', num2str(win_size(2)),'win_',num2str(smooth_size),'smooth'];

    fig_par_name = [num2str(smooth_size), 'smooth_', num2str(win_size(1)), 'to', num2str(win_size(2)), 'win'];
    
else
    
    display('win_size must be a scalar (lower limit of beta epoch length) or an interval (lower and upper limits).')
    
    return
    
end

% record_label = {'st_m1', 'st_stn', 'st_m1_6OHDA'}; record_periods = {1:2, 1:2, 1}; record_channels = [2, 1, 2];
% 
% fig_legend = {'Motor Beta, Pre-Infusion', 'Motor Beta, Post-Infusion', 'STN Beta, Pre-Infusion', 'STN Beta, Post-Infusion', 'Motor Beta, 6OHDA'};
% 
% color_map = [0 0 1; 0 .5 0; 1 0 1; 0 1 1; 1 0 0];

record_label = {'st_m1', 'st_m1_6OHDA'}; record_periods = {1:2, 1}; record_channels = [2, 2];

fig_legend = {'Motor Beta, Pre-Infusion', 'Motor Beta, Post-Infusion', 'Motor Beta, 6OHDA'};

color_map = [0 0 1; 0 .5 0; 1 0 0];

no_records = length(record_label);

figure;

[All_beta_length_histograms, All_beta_length_bins] = deal(nan(no_bins, 2*no_records));

max_beta_length = 0;

min_beta_length = 3*10^6;

Beta_lengths = cell(2, no_records);

for r = 1:no_records

    %% Collecting beta epoch lengths, pre and post.
    
    load([record_label{r}, '_subjects.mat'])
    
    load([record_label{r}, '_beta_', par_name, '_length_hist.mat'], 'All_beta_lengths')
    
    for pd = 1:length(record_periods{r})
    
        Beta_lengths{pd, r} = All_beta_lengths{record_periods{r}(pd), record_channels(r)};
        
        max_beta_length = max(max_beta_length, max(Beta_lengths{pd, r}));
        
        min_beta_length = min(min_beta_length, min(Beta_lengths{pd, r}));
        
    end
    
end

bin_edges = exp(linspace(log(min_beta_length), log(max_beta_length + 1), no_bins + 1));

bin_centers = sqrt(bin_edges(1 : end - 1).*bin_edges(2 : end));

for r = 1:no_records
    
    for pd = 1:length(record_periods{r})
        
        if max_beta_length > 0
            
            [h, ~] = histc(Beta_lengths{record_periods{r}(pd), r}, bin_edges);
            
            All_beta_length_histograms(:, (r - 1)*2 + pd) = h(1 : end - 1);
            
            All_beta_length_bins(:, (r - 1)*2 + pd) = bin_centers;
            
        end
        
    end
    
end

All_beta_length_histograms(All_beta_length_histograms == 0) = eps;

set(0, 'DefaultAxesColorOrder', color_map)

loglog(All_beta_length_bins, All_beta_length_histograms, '-o')

axis tight

ylim([1 - eps, max(max(All_beta_length_histograms))])

legend(fig_legend)

title(['Distribution of \beta Epochs (\beta > 2 S.D.) by Length, for ', num2str(smooth_size/1000), ' s Smoothing'])

xlabel('Epoch Length')

ylabel('Number Observed')

save_as_pdf(gcf, ['Carb_6OHDA_', fig_par_name, 'beta_length_hist'])

end