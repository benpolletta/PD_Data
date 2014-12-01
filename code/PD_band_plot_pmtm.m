function PD_band_plot_pmtm(subjects_mat, epoch_secs, window_secs)

close('all')

load(subjects_mat)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[r, c] = subplot_size(no_bands);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_BP_windowed.mat'])
    
    t_win_BP = t_win;
    
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_stats.mat'], 't_win', 'p_less', 'p_greater')
        
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.mat'])
    
    lower_test = p_less <= .01/(2*sum_all_dimensions(~isnan(p_less))); lower_test = +lower_test;
    
    lower_test(lower_test == 0) = nan;
    
    upper_test = p_greater <= .05/(2*sum_all_dimensions(~isnan(p_greater))); upper_test = +upper_test;
    
    upper_test(upper_test == 0) = nan;
    
    t = All_data.t;
    
    for n = 1:no_norms
        
        %% Plotting Band Power.
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        BP_data(artifact_indicator, :, :) = nan;
        
        figure((fo - 1)*no_norms + n)
        
        for b = 1:no_bands
            
            subplot(r, c, b);
            
            [ax, ~, ~] = plotyy(t, reshape(BP_data(:, b, :), size(BP_data, 1), 2), t_win_BP, median_pow(:, :, b, n));
            
            hold(ax(1), 'on'), hold(ax(2), 'on')
            
            axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
            
            set(ax,{'ycolor'},{'black';'black'})
            
            set(get(ax(1),'YLabel'), 'String', ['Power', norm_labels{n}])
            
            set(get(ax(2),'YLabel'), 'String', 'Mean Power')
            
            add_stars(ax(1), t_win, [lower_test(:, :, b, n) upper_test(:, :, b, n)], [0 0 1 1], [0 0 1; 0 .5 0; 0 0 1; 0 .5 0])
            
            xlabel('Time (s)')
            
            ylabel(['Power', norm_labels{n}])
            
            title([folder, ', ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
            
            if b == 1
            
                legend(ax(1), chan_labels, 'Location', 'NorthWest')
            
            end
            
        end
        
        try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_' , num2str(window_secs), 's_stats_pmtm_BP', norms{n}]), end
        
    end
    
    % close('all')
    
end