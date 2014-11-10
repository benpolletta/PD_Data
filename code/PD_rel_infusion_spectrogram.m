function PD_rel_infusion_spectrogram(subject_mat, chunk_size)

load(subject_mat)

sampling_freq = 1000;

freqs = 1:200; no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 100 200; 0 200]; no_bands = size(bands, 1);

for b = 1:no_bands, band_indices = freqs <= bands(b, 2) & freqs >= bands(b, 1); end

no_cycles = linspace(3, 21, no_freqs);

beta_indices = freqs >= 8 & freqs <= 30;

no_channels = length(chan_labels);

suffixes = {'', '_norm', '_BP', '_BP_norm'}; no_suffixes = length(suffixes);

for s = 1:no_suffixes
   
    if s <= no_suffixes/2
    
        formats{s} = make_format(no_freqs, 'f');
        
    else
        
        formats{s} = make_format(no_bands, 'f');
    
    end
    
end

BP_chunk = nan(chunk_size*sampling_freq, no_bands);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    data_length = length(PD_dec);
    
    no_mins = min(floor(data_length/(sampling_freq*60)), (basetime + 1500)/60);
    
    no_chunks = no_mins*60/chunk_size;
    
    for s = 1:no_suffixes
        
        for ch = 1:no_channels
       
            fid_mat(s, ch) = fopen([subj_name, '_ch', num2str(ch), '_wt_', suffix{s}, '.txt'], 'w');
        
        end
        
    end
    
    for c = 1:no_chunks
        
        chunk_name = [subj_name, '_chunk', num2str(c), '_wt'];
    
        chunk_start = (c - 1)*chunk_size*sampling_freq + 1;
        
        chunk_end = c*chunk_size*sampling_freq;
        
        t = ((chunk_start:chunk_end)/sampling_freq - basetime);
        
        minute = ceil(mean(t)/60);
        
        figure
        
        for ch = 1:no_channels
            
            LFP_chunk = PD_dec(chunk_start:chunk_end, ch);
            
            %% Calculate Spectrogram and Band Power.
            
            WT_chunk = abs(wavelet_spectrogram(LFP_chunk, sampling_freq, freqs, no_cycles));
            
            for b = 1:no_bands
            
                BP_chunk(:, b) = sqrt(sum(WT_chunk(:, band_indices{b}, ch).^2, 2));
            
            end
            
            BP_norm = BP_chunk./repmat(BP_chunk(:, end), no_bands);
            
            WT_norm = WT_chunk./repmat(BP_chunk(:, end), no_freqs);
            
            %% Print Spectrogram and Band Power.
            
            fprintf(fid_mat(1, ch), formats{1}, WT_chunk');
            
            fprintf(fid_mat(2, ch), formats{2}, WT_norm');
            
            fprintf(fid_mat(3, ch), formats{3}, BP_chunk');
            
            fprintf(fid_mat(4, ch), formats{4}, BP_norm');
            
            %% Plot Spectrogram.
            
            subplot(5, 1, (ch - 1)*2 + 1)
            
            imagesc(t, freqs, WT_chunk(:, :, ch)')
            
            if ch == 1
               
                title(['Wavelet Spectrogram, ', folder, ', Minute ', num2str(minute)])
                
            end
            
            ylabel(chan_labels{ch})
            
            axis xy
            
            hold on, plot(t', ones(length(t), 2)*diag([10 30]), ':w'),
            
            subplot(5, 1, (ch - 1)*2 + 2)
            
            imagesc(t, freqs, zscore(WT_chunk(:, :, ch))')
            
            ylabel('Freq. Normalized')
            
            axis xy
            
            hold on, plot(t', ones(length(t), 2)*diag([10 30]), '--w'),
            
        end
        
        %% Plot LFP.
        
        subplot(5, 1, 5)
        
        [ax, h1, h2] = plotyy(t, LFP_chunk, t, beta_chunk);
        
        %title(['Raw LFP & Pct. \beta Power, ', folder, ', Minute ', num2str(minute)])
        
        %box off
        
        axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1))) %axis tight
        
        set(ax,{'ycolor'},{'black';'black'})
        
        set(get(ax(1),'YLabel'),'String','LFP')
        
        set(get(ax(2),'YLabel'),'String','% \beta Power')
        
        legend(h1, chan_labels, 'Location', 'NorthWest')
        
        legend(h2, chan_labels, 'Location', 'SouthWest')
        
        save_as_pdf(gcf, chunk_name)
        
    end
    
    close('all')
    
    % end
    
end

end