function PD_rel_infusion_plot_spectrogram(subject_mat, chunk_size)

load(subject_mat)

sampling_freq = 1000;

freqs = 1:200; no_freqs = length(freqs);

no_cycles = linspace(3, 21, no_freqs);

bands = [1 8; 8 30; 30 100; 100 200; 0 200]; no_bands = size(bands, 1);

[band_indices, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands

    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    band_labels{b} = [num2str(bands(b, 1)), '-', num2str(bands(b, 2)),'Hz'];

end
    
no_channels = length(chan_labels);

BP_chunk = nan(chunk_size*sampling_freq, no_bands); %[total_chunk, beta_chunk] = deal(nan(chunk_size*sampling_freq, 2));

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    PD_zs = zscore(PD_dec);
    
    data_length = length(PD_dec);
    
    no_mins = min(floor(data_length/(sampling_freq*60)), (basetime + 1500)/60);
    
    no_chunks = no_mins*60/chunk_size;
    
    load([subj_name,'_wt.mat'], 'Spec', 'BP')
    
    BP_zs = nan(size(BP));
    
    for ch = 1:2
    
        BP_zs(:, :, ch) = zscore(BP(:, :, ch));
    
    end
        
    for c = (19*3):no_chunks
        
        chunk_name = [subj_name, '_chunk', num2str(c), '_wt'];
    
        chunk_start = (c - 1)*chunk_size*sampling_freq + 1;
        
        chunk_end = c*chunk_size*sampling_freq;
        
        LFP_chunk = PD_zs(chunk_start:chunk_end, :);
        
        WT_chunk = abs(Spec(chunk_start:chunk_end, :, :)); % WT_chunk = abs(wavelet_spectrogram(LFP_chunk, sampling_freq, freqs, no_cycles));
        
        t = ((chunk_start:chunk_end)/sampling_freq - basetime);
        
        minute = ceil(mean(t)/60);
        
        figure
        
        for ch = 1:no_channels
            
            % %% Compute Band Power.
            % 
            % for b = 1:no_bands
            % 
            %     BP_chunk(:, b) = sqrt(sum(WT_chunk(:, band_indices{b}, ch).^2, 2));
            % 
            % end
            
            %% Plot Spectrogram.
            
            subplot(6, 1, (ch - 1)*3 + 1)
            
            imagesc(t, freqs, WT_chunk(:, 1:100, ch)')
            
            if ch == 1
               
                title(['Wavelet Spectrogram, ', folder, ', Minute ', num2str(minute)])
                
            end
            
            ylabel(chan_labels{ch})
            
            axis xy
            
            hold on, plot(t', ones(length(t), 2)*diag([10 30]), ':w'),
            
            subplot(6, 1, (ch - 1)*3 + 2)
            
            imagesc(t, freqs, zscore(WT_chunk(:, 1:100, ch))')
            
            ylabel('Freq. Normalized')
            
            axis xy
            
            hold on, plot(t', ones(length(t), 2)*diag([10 30]), '--w'),
            
            %% Plot LFP & Band Power.
            
            subplot(6, 1, (ch - 1)*3 + 3)
            
            plot(t, [LFP_chunk(:, ch) BP_zs(chunk_start:chunk_end, :, ch)])%BP_chunk]))

            axis tight
            
            yl = ylim;
            
            set(gca, 'YTick', ceil(yl(1)):floor(yl(2)), 'YTickLabel', ceil(yl(1)):2:floor(yl(2)), 'YGrid', 'on')
            
            ylabel(['LFP & Band Power'])
            
            if ch == 1
            
                legend({'LFP', band_labels{:}})
            
            end
            
        end
        
        % %% Plot LFP.
        % 
        % subplot(5, 1, 5)
        % 
        % [ax, h1, h2] = plotyy(t, LFP_chunk, t, beta_chunk./total_chunk);
        % 
        % %title(['Raw LFP & Pct. \beta Power, ', folder, ', Minute ', num2str(minute)])
        % 
        % box off
        % 
        % axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1))) %axis tight
        % 
        % set(ax,{'ycolor'},{'black';'black'})
        % 
        % set(get(ax(1),'YLabel'),'String','LFP')
        % 
        % set(get(ax(2),'YLabel'),'String','% \beta Power')
        % 
        % legend(h1, chan_labels, 'Location', 'NorthWest')
        % 
        % legend(h2, chan_labels, 'Location', 'SouthWest')
        
        save_as_pdf(gcf, chunk_name)
        
    end
    
    close('all')
    
    % end
    
end

end