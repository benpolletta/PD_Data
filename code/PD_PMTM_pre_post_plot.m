function PD_PMTM_pre_post_plot(filename, folder_indices, band_index, significance)

figure

pd_names = {'Pre-Infusion', 'Post-Infusion'}; no_periods = length(pd_names);

% pd_linestyles = {'--', '-'};

channel_labels = {'Striatum'; 'Motor Cortex'};

no_channels = length(channel_labels);

analysis_name = make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'pmtm',{[150 150],[1 1]},2);

load(analysis_name)

sampling_freq = 500;

f = sampling_freq*(1:size(SW, 1))/(2*size(SW,1));

if isempty(folder_indices)
    
    folder_indices = ones(size(SW, 3), 1);
    
    folder_indices([3 9]) = 0; folder_indices = logical(folder_indices);
    
    folder_indices = {'', folder_indices};
    
end

stats_name = [make_sliding_window_analysis_name([filename, folder_indices{1},...
    '_pre_post_band', num2str(band_index)], 'pmtm', {[150 150],[1 1]}, 2), '_ttest'];

load(stats_name)

test = p_vals < significance;

% SW_mean(:, :, :, :) = squeeze(nanmean(SW(:, :, folder_indices{2}, :), 3));
% SW_se(:, :, :, :) = squeeze(nanstd(SW(:, :, folder_indices{2}, :), [], 3)/sqrt(sum(folder_indices{2})));

SW_post_by_pre = 100*SW(:, :, folder_indices{2}, 2)./SW(:, :, folder_indices{2}, 1) - 100;

SW_mean = squeeze(nanmean(SW_post_by_pre, 3));
SW_se = squeeze(nanstd(SW_post_by_pre, [], 3)/sqrt(sum(folder_indices{2})));

SW_pre_from_post = squeeze(diff(SW(:, :, folder_indices{2}, :), [], 4));

freq_limit = 50;
freq_indices = f <= 50;

for channel = 1:no_channels
    
    h = subplot(2, no_channels, channel);
    
    boundedline(f(freq_indices), squeeze(SW_mean(freq_indices, channel)),...
        prep_for_boundedline(squeeze(SW_se(freq_indices, channel)))); 
        % squeeze(nanmean(SW_pre_from_post(freq_indices, :, channel), 2)),... 
        % prep_for_boundedline(squeeze(nanstd(SW_pre_from_post(freq_indices, :, channel), [], 2)))); % 
    
    hold on
            
    add_stars(h, f(freq_indices), squeeze(test(freq_indices, channel, :)), [1 0], [1 0 0; 1 .5 0])
    
    set(gca, 'FontSize', 14)
    
    title(channel_labels{channel})
    
    legend(pd_names)
    
    ylabel({'Spectral Power'; 'Mean \pm S.E.'})
    
    h = subplot(2, no_channels, no_channels + channel);
        
    set(h, 'NextPlot', 'add', 'ColorOrder', distinguishable_colors(sum(folder_indices{2})), 'FontSize', 14)
    
    % for pd = 1:2
    
        plot(f(freq_indices), nanzscore(squeeze(SW_post_by_pre(freq_indices, channel, :))')') % , pd_linestyles{pd})
        
    % end
    
    hold on
    
    ylabel({'Spectral Power'; 'Post-Infusion - Pre-Infusion'})
    
    xlabel('Freq. (Hz)')
    
end

fig_name = make_sliding_window_analysis_name([filename, folder_indices{1},...
    '_pre_post_band', num2str(band_index), '_p', num2str(significance)],...
    'pmtm',{[150 150],[1 1]},2);

save_as_pdf(gcf, fig_name)

