function PD_boundedline_collect_pmtm(subjects_mat, epoch_secs, window_secs)

close('all')

load(subjects_mat), no_subjects = length(folders);

load([folders{1}, '/', prefixes{1}, '_', num2str(epoch_secs), 's_epoch_pmtm.mat'], 'freqs')

load([folders{1}, '/', prefixes{1}, '_all_channel_data_dec.mat'], 'sampling_freq')

plot_freq_indices = freqs >= 0 & freqs <= 200; plot_freqs = freqs(plot_freq_indices);

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

window_length = window_secs*sampling_freq;

epochs_per_window = window_secs/epoch_secs;

All_pow = nan(no_subjects*epochs_per_window, length(plot_freqs), 8, 2, no_bands, no_norms);

[All_mean, All_std] = deal(nan(no_subjects, length(plot_freqs), 8, 2, no_bands, no_norms));

All_indices = nan(4, 2, no_subjects, no_bands, no_norms);

for fo = 1:no_subjects
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_BP_windowed.mat'])
    
    t_win_BP = t_win;
    
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_stats.mat'], 't_win', 'p_less', 'p_greater')
        
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.mat'])
    
    t = All_data.t;
    
    for b = 1:no_bands
        
        for n = 1:no_norms
            
            %% Collecting Spectra.
            
            Spec_data = getfield(All_data, ['Spec', norms{n}]);
            
            Spec_data(artifact_indicator, :, :) = nan;
                
            indices = get_indices_pow_p(t_win_BP, median_pow(:, :, b, n), p_less(:, :, b, n), p_greater(:, :, b, n));
                
            All_indices(:, :, fo, n, b) = indices;
            
            for ch = 1:2
                
                Spec_window = nan(epochs_per_window, length(plot_freqs), 8);
                
                [mean_window, std_window] = deal(nan(length(plot_freqs), 8));
                
                for i = 1:size(indices, 1)
                    
                    win_indices = t >= (indices(i, ch) - 1)*epoch_secs & t < (indices(i, ch) + epochs_per_window)*epoch_secs;
                    
                    Sw = Spec_data(win_indices, plot_freq_indices, :);
                    
                    % size(Sw)
                    
                    Spec_window(1:min(24, size(Sw, 1)), :, (2*i - 1):(2*i)) = Sw(1:min(24, size(Sw, 1)), :, :);
                    
                    mean_window(:, (2*i - 1):(2*i)) = reshape(nanmean(Sw), length(plot_freqs), 2);
                    
                    std_window(:, (2*i - 1):(2*i)) = reshape(nanstd(Sw), length(plot_freqs), 2);
                    
                end
                
                All_pow((fo - 1)*epochs_per_window + (1:epochs_per_window), :, :, ch, b, n) = Spec_window;
                
                All_mean(fo, :, :, ch, b, n) = nanmean(Spec_window);
                
                All_std(fo, :, :, ch, b, n) = nanstd(Spec_window);
                
            end
            
        end
        
    end
    
end

save([subjects_mat(1:(end - length('subjects.mat'))), '_boundedline_collect_pmtm.mat'], 'All_indices', 'All_pow', 'All_mean', 'All_std')

end

function indices = get_indices_pow_p(t_win_pow, power, p_less, p_greater)

[~, increase_pre] = nanmax(power(t_win_pow < 0, :));

[~, decrease_pre] = nanmin(power(t_win_pow < 0, :));

[~, increase_post] = nanmax(p_greater(:, :));

increase_post = increase_post + sum(t_win_pow <= 0);

[~, decrease_post] = nanmax(p_less(:, :));

decrease_post = decrease_post + sum(t_win_pow <= 0);

indices = [increase_pre; decrease_pre; increase_post; decrease_post];

end