function PD_stats_pmtm(subjects_mat, epoch_secs, window_secs)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

load(subjects_mat)

epochs_per_window = window_secs/epoch_secs;

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    t = All_data.t;
    
    no_windows = floor((t(end) - window_secs)/epoch_secs);
        
    t_win = nan(no_windows, 1);
    
    [p_less, p_greater] = deal(nan(no_windows, 2, no_bands, no_norms));
    
    for n = 1:no_norms
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        %% Eliminating segments with high total power.
        
                    
        %% Comparing Band Power against Baseline, for Each Window.

        for b = 1:no_bands
            
            for w = 1:no_windows
                
                win_indices = t >= (w - 1)*epoch_secs & t < (w + epochs_per_window)*epoch_secs;
                
                t_win(w) = median(t(win_indices));
                
                for ch = 1:2
                    
                    p_less(w, ch, b, n) = ranksum(BP_data(t < 0, b, ch), BP_data(win_indices, b, ch), 'tail', 'right');
                    
                    p_greater(w, ch, b, n) = ranksum(BP_data(t < 0, b, ch), BP_data(win_indices, b, ch), 'tail', 'left');
                    
                end
                
            end
            
        end
        
    end
    
    save([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_stats.mat'], 't_win', 'bands', 'norms', 'p_less', 'p_greater')
    
end