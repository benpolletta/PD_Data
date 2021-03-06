function PD_GC_pre_post_figure_new(significance)

filename = 'STR_w_M1';

band_index = 4;

if isempty(significance), significance = .025; end

figure

% pd_names = {'Post-Infusion - Pre-Infusion'}; % , 'Post-Infusion'}; 

% no_periods = length(pd_names);

% pd_linestyles = {'--', '-'};

direction_indices = {{2, 1, ':'}, {1, 2, ':'}};

direction_labels = {'\bf{Striatum \rightarrow Motor Cortex}'; '\bf{Motor Cortex \rightarrow Striatum}'};

no_directions = length(direction_labels);

varargin = {[], '', 1};

analysis_name = make_sliding_window_analysis_name([filename, '_pre_post_band', num2str(band_index)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2, varargin{:});

load(analysis_name)

GC = SW;

sampling_freq = 500;

f = sampling_freq*(1:size(GC, 1))/(2*size(GC,1));

freq_limit = 50;

freq_indices = find(f <= 50);

direction_indices = {{freq_indices, 1, 2, ':'}, {freq_indices, 2, 1, ':'}};

folder_index_cell = make_folder_index_cell;

no_folder_indices = length(folder_index_cell);

group_titles = {{'{\bf All Data}';'{\fontsize{10}Partial Directed Coherence}';...
    '{\fontsize{10}Post-Infusion - Pre-Infusion}';'{\fontsize{10}Mean \pm SE}'},...
    {'{\bf M1+}';'{\fontsize{10}\Delta PDC (Mean \pm SE)}'}, {'{\bf M1-}';'{\fontsize{10}\Delta PDC (Mean \pm SE)}'}};

for fi = 1:no_folder_indices
    
    folder_indices = folder_index_cell{fi};
    
    % PD_GC_pre_post_ttest(filename, folder_indices, band_index)
    
    stats_name = [make_sliding_window_analysis_name([filename, folder_indices{1},...
        '_pre_post_band', num2str(band_index)], 'mvgc_analysis', {[150 150], [2 2]},...
        2, varargin{:}), '_ttest'];
    
    load(stats_name)
    
    % GC_post_over_pre = squeeze(100*GC(:, :, :, folder_indices{2}, 2)./GC(:, :, :, folder_indices{2}, 1) - 100);
    
    % GC_post_over_pre(:, :, :, :, 2) = 0;
    
    GC_pre_from_post = squeeze(diff(GC(:, :, :, folder_indices{2}, :), [], 5));
    
    GC_pre_from_post(:, :, :, :, 2) = 0;
    
    GC_plot = GC_pre_from_post; % GC_plot = GC_post_over_pre; % GC_plot = GC(:, :, :, folder_indices{2}, :); % 
    
    GC_mean = squeeze(nanmean(GC_plot, 4));
    GC_se = squeeze(nanstd(GC_plot, [], 4)/sqrt(sum(folder_indices{2})));
    GC_ci = norminv(1 - significance)*GC_se;
    
    for direction = 1:no_directions
        
        ax(fi, direction) = subplot(no_folder_indices, no_directions, (fi - 1)*no_directions + direction);
    
        % set(ax(fi, direction), 'NextPlot', 'add', 'ColorOrder', [0 0 0; 0 0 0])
        
        boundedline(f(freq_indices), diag((f(freq_indices)).^(2/3))*squeeze(GC_mean(direction_indices{direction}{:})),...
            prep_for_boundedline(diag((f(freq_indices)).^(2/3))*squeeze(GC_ci(direction_indices{direction}{:}))));
        % squeeze(nanmean(GC_pre_from_post(freq_indices, :, direction), 2)),...
        % prep_for_boundedline(squeeze(nanstd(GC_pre_from_post(freq_indices, :, direction), [], 2)))); %
        
        axis tight
        
        % ylim([0 .75]) % ylim([-1 1]*1000)
        
        hold on
        
    end

    sync_axes(ax(fi, :)) % linkaxes(fliplr(ax(fi, :))), linkaxes(ax(fi, :), 'off')
    
    for direction = 1:no_directions
        
        subplot(no_folder_indices, no_directions, (fi - 1)*no_directions + direction)
        
        % plot(f(freq_indices), zeros(size(f(freq_indices))), 'k')
    
        test = squeeze(p_vals(direction_indices{direction}{:}) < significance);
        
        add_stars(ax(fi, direction), f(freq_indices), test, [1 0], [1 0 0; 1 .5 0]) % test(freq_indices, :, direction), [1 0], [1 0 0; 1 .5 0])
        
        xlim([0 freq_limit])
        
        set(gca, 'FontSize', 14)
        
        if fi == 1
            
            title(direction_labels{direction})
            
        end
        
        % if fi == 1 && direction == 1
        % 
        %     legend(pd_names)
        % 
        % end
        
        if direction == 1
        
            ylabel(group_titles{fi})
            
        end
        
        if fi == no_folder_indices
        
            xlabel('Freq. (Hz)')
            
        end
    
    end
    
end

fig_name = make_sliding_window_analysis_name([filename,...
    '_pre_post_band', num2str(band_index), '_p', num2str(significance)],...
    'mvgc_analysis', {[150 150], [2 2]}, 2, varargin{:});

save_as_pdf(gcf, [fig_name, '_altogether'])

end

function folder_index_cell = make_folder_index_cell

folder_indices = ones(14, 1);

folder_indices([3 9]) = 0; folder_indices = logical(folder_indices);

folder_indices = {'', folder_indices};

folder_index_cell{1} = folder_indices;

load('M1_groups')

folder_index_cell{2} = M1_increased_index;

folder_index_cell{3} = M1_not_increased_index;

end

