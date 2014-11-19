function PD_stats_pmtm(subjects_mat, epoch_length)

load(subjects_mat)

sampling_freq = 1000;

%freqs = 1:200;

bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]; no_bands = size(bands, 1);

[r, c] = subplot_size(no_bands);

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_length/sampling_freq), 's_epoch_pmtm.mat']);
    
    t = All_data.t;
    
    freqs = All_data.freqs;
    
    plot_freq_indices = freqs >= 0 & freqs <= 200;
    
    for n = 1:no_norms
        
        %% Comparing Band Power against Baseline.
        
        BP_data = getfield(All_data, ['BP', norms{n}]);
        
        figure((fo - 1)*no_norms*2 + no_norms + n)
        
        
            
            plot(t, zscore(reshape(BP_data(:, b, :), size(BP_data, 1), 2)))
            
            axis tight
            
            xlabel('Time (s)')
            
            ylabel(['Power', norm_labels{n}])
            
            title(['Power, ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
            
            legend(chan_labels)
            
        end
        
        try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_length/sampling_freq), 's_pmtm_BP', norms{n}]), end
        
    end
    
    % close('all')
    
end