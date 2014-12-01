function PD_colorplot_pmtm(subjects_mat, epoch_secs)

close('all')

load(subjects_mat)

norms = {'', '_norm', '_pct', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Total Power', ', % Baseline Power', ', Baseline Normalized % Total Power'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    t = All_data.t;
    
    freqs = All_data.freqs;
    
    plot_freq_indices = freqs >= 0 & freqs <= 200;
    
    for n = 1:no_norms
            
        %% Plotting Spectrograms.
        
        Spec_data = getfield(All_data, ['Spec', norms{n}]);
        
        Spec_data(artifact_indicator, :, :) = nan;
        
        for ch = 1:2
            
            figure((fo - 1)*no_norms*2 + n)
            
            subplot(2, 1, ch)
            
            if n < 3
                
                Spec_plot = nanzscore(Spec_data(:, plot_freq_indices, ch))';
                
            else
                
                Spec_plot = Spec_data(:, plot_freq_indices, ch)';
                
            end
            
            imagesc(t, freqs(plot_freq_indices), Spec_plot)
            
            cl = caxis; caxis([cl(1) median(max(Spec_plot, [], 2))])
            
            axis xy
            
            xlabel('Time (s)'); ylabel('Hz');
            
            title(['Multi-taper Spectrogram of ', folder, ', ', chan_labels{ch}, norm_labels{n}])
            
            grid on
            
        end
        
        try, save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_epochs_pmtm', norms{n}]), end
        
    end
    
    % close('all')
    
end