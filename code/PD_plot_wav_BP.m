function PD_plot_wav_BP(subject_mat)

% if isscalar(win_size)
% 
%     par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];
% 
% elseif length(win_size) == 2
%    
%     par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size(1)), 'to', num2str(win_size(2)),'win_',num2str(smooth_size),'smooth'];
%     
% else
%     
%     display('win_size must be a scalar (lower limit of beta epoch length) or an interval (lower and upper limits).')
%     
%     return
%     
% end
    
load(subject_mat)

sampling_freq = 1000;

norms = {'', '_pct', '_norm', '_norm_pct'}; no_norms = length(norms);

norm_labels = {'', ', % Baseline', ', % Total Power', ', % Baseline Total Power'};

for n = 1:3 %no_norms
    
    for fo = 1:length(folders)
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        basetime = basetimes(fo);
        
        % base_index = basetime*sampling_freq;
        % 
        % no_pre_epochs = floor(base_index/epoch_length);
        % 
        % start_index = base_index - no_pre_epochs*epoch_length;
        % 
        % no_post_epochs = floor((size(PD_dec, 1) - basetime*sampling_freq)/epoch_length);
        % 
        % no_epochs = no_pre_epochs + no_post_epochs;
        
        BP_data = load([subj_name,'_wt_BP.mat'], ['BP', norms{n}], 'bands');
        
        bands = getfield(BP_data, 'bands'); no_bands = length(bands); [r, c] = subplot_size(no_bands);
        
        for b = 1:no_bands, band_labels{b} = sprintf('%.0f - %.0f Hz', bands(b, :)); end
        
        BP_data = getfield(BP_data, ['BP', norms{n}]);
        
        t = (1:size(BP_data, 1))/sampling_freq - basetime;
        
        %% Plotting Band Power.
        
        figure((2*n - 2)*length(folders) + fo)
        
        max_BP = max(max(max(BP_data))); min_BP = min(min(min(BP_data)));
        
        for b = 1:no_bands
        
            subplot(r, c, b)
        
            plot(t, reshape(BP_data(:, b, :), size(BP_data, 1), size(BP_data, 3)))
        
            axis tight
            
            ylim([min_BP, max_BP])
        
            xlabel('Time (s)')
        
            ylabel(['Power', norm_labels{n}])
        
            title([folder, ', Power, ', num2str(bands(b, 1)), ' - ', num2str(bands(b, 2)), ' Hz'])
        
            if b == 1
        
                legend(chan_labels)
        
            end
        
        end
        
        try save_as_pdf(gcf, [subj_name, '_wav_BP', norms{n}, '_subplots']), end
        
        BP_data = zscore(BP_data);
        
        figure((2*n - 1)*length(folders) + fo) %(fo - 1)*2*no_norms + 2*n)
        
        for ch = 1:2
           
            subplot(2, 1, ch)
            
            plot(t, BP_data(:, :, ch))
            
            axis tight
            
            if ch == 1
                
                title([folder, ', Power', norm_labels{n}])
                
                legend(band_labels, 'Location', 'NorthWest')
               
            else
                
                xlabel('Time (s)')
                
            end
            
            ylabel(chan_labels{ch})
            
        end
        
        try save_as_pdf(gcf, [subj_name, '_wav_BP', norms{n}]), end
            
        % t_mat = repmat(t, 1, size(BP_data, 2));
        % 
        % hold on
        % 
        % plot(t_mat(logical(BP_high(:, :, ch))), BP_data(logical(BP_high(:, :, ch)), :, ch), 'LineWidth', 2)
        % 
        % plot(t(BP_high_cum(:, :, ch)), BP_data(BP_high_cum(:, :, ch), :, ch), ':', 'LineWidth', 2)
        % 
        % plot(t(BP_high(:, :, ch) & BP_high_cum(:, :, ch)), BP_data(BP_high(:, :, ch) & BP_high_cum(:, :, ch), :, ch), 'LineWidth', 3)
        
    end
          
end