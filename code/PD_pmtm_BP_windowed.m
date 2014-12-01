function PD_pmtm_BP_windowed(subjects_mat, epoch_secs, window_secs)

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

load(subjects_mat)

epochs_per_window = window_secs/epoch_secs;

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
        
    load([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.mat'])
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    t = All_data.t;
    
    t(artifact_indicator) = nan;
    
    no_windows_pre = floor((abs(nanmin(t(t < 0))) - window_secs)/epoch_secs);
    
    start_sec = -(no_windows_pre*epoch_secs + window_secs);
    
    no_windows_post = floor((nanmax(t) - window_secs)/epoch_secs);
    
    no_windows = no_windows_pre + no_windows_post;
    
    t_win = nan(no_windows, 1);
    
    [median_pow, mean_pow, std_pow] = deal(nan(no_windows, 2, no_bands, no_norms));
    
    for n = 1:no_norms
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        BP_data(artifact_indicator, :, :) = nan;
                    
        %% Comparing Band Power against Baseline, for Each Window.

        for b = 1:no_bands
            
            for w = 1:no_windows
                
                win_start = (w - 1)*epoch_secs + start_sec;
                
                win_end = (w + epochs_per_window)*epoch_secs + start_sec;
                
                win_indices = t >= win_start & t < win_end;
                
                t_win(w) = median(t(win_indices));
                
                for ch = 1:2
                    
                    median_pow(w, ch, b, n) = nanmedian(BP_data(win_indices, b, ch));
                    
                    mean_pow(w, ch, b, n) = nanmean(BP_data(win_indices, b, ch));
                    
                    std_pow(w, ch, b, n) = nanstd(BP_data(win_indices, b, ch));
                    
                end
                
            end
            
        end
        
    end
    
    save([subj_name, '_', num2str(epoch_secs), 's_pmtm_', num2str(window_secs), 's_BP_windowed.mat'], 't_win', 'bands', 'norms', 'median_pow', 'mean_pow', 'std_pow')
    
end