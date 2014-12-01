function [mean_E, std_E, normalizer] = PD_powerfreq_v2(subject_mat)

low_freq_lim = 3;
high_freq_lim = 72;
no_freqs = high_freq_lim - low_freq_lim + 1;

load(subject_mat)
no_folders = length(folders);
no_channels = length(chan_labels);

pd_label = {'pre','post'};

[normalizer, avg_E] = deal(nan(no_freqs, no_folders, no_channels));

[mean_E, std_E] = deal(nan(no_freqs, no_folders, 2));

for pd=1:length(pd_label)
    
    clear avg_E pre_E
    
    for fo = 1:no_folders % - 1
        
        folder = folders{fo};
        prefix = prefixes {fo};
        
        load([folder, '/', prefix, '_all_channel_data_dec.mat'], 'sampling_freq')
        
        outputname=[subject_mat(1:end-12), prefix, '_', pd_label{pd}, '_all_channel_spec_HT', '.mat'];
        
        load(outputname)
        
        for ch = 1:no_channels
            
            % pre_E = zeros(no_freqs, length(energy), no_channels);
            % 
            % pre_E(:, :, ch) = energy(low_freq_lim:high_freq_lim, :, ch);
            
            % %% Using beta epochs to figure out when to start looking at spectrogram.
            %
            % beta_listname = [folder,'/',prefix,'_ch1_beta_post_7out_2sd_333.3333win_5000smooth_win.list'];
            % 
            % betapower = [folder,'/',prefix,'_all_channel_data_dec_HAP.mat'];
            % 
            % % Loading block numbers, epoch start and end indices.
            % [blocks, beta_starts, beta_ends, epochs, epoch_starts, epoch_ends] = text_read(beta_listname,'%f%f%f%f%f%f%*[^\n]');
            % 
            % % find(max(epochs)
            % betatime = beta_starts(find(epochs == max(epochs)) - max(epochs) + 1);
            
            % pre_E(:, fo, ch) = mean(zscore(log(energy(a:b, :, ch)), 0, 1), 2); % Taking mean zscore of log of energy at each frequency.
            
            if pd==1 % Using first 30 seconds of pre-infusion data to normalize.
                
                normalizer(:, fo, ch) = mean(energy(low_freq_lim:high_freq_lim, 1:30*sampling_freq, ch), 2);
                % normalizer(:, fo, ch) = mean(energy(low_freq_lim:high_freq_lim, :, ch), 2);
                
            end
            
            pre_E = zeros(no_freqs, length(energy), no_channels);
            
            % pre_E(:,:,ch) = energy(a:b, :, ch) ./ repmat(normalizer(:, fo, ch), [1 length(energy) 1]); % Divisive normalization.
            
            pre_E(:, :, ch) = (energy(low_freq_lim:high_freq_lim, :, ch) - repmat(normalizer(:, fo, ch), [1 length(energy) 1])) ./ ...
                (energy(low_freq_lim:high_freq_lim, :, ch) + repmat(normalizer(:, fo, ch), [1 length(energy) 1])); % Normalization as a kind of modulation index.
            
            % avg_E(:,ch) = reshape(pre_E(:, fo, ch), b - a + 1, 1); % Renaming pre_E from above.
            
            figure(fo)
            
            subplot(1,2,ch)
            
            if pd == 1
                
                avg_E(:, fo, ch) = mean(pre_E(:, 1:30*sampling_freq, ch), 2); % Averaging normalized data over first 30 seconds.
                
                plot(avg_E(:, fo, ch)) % Plotting average.
                
                hold on
                
                figure(fo + 10)
                
                subplot(2, 2, ch) % Plotting all 30 seconds (normalized).
                
                imagesc(1:30*sampling_freq, 1:size(pre_E, 1), pre_E(:, 1:30*sampling_freq, ch))
                
                axis xy
            
            else
                
                % avg_E(:, fo, ch) = mean(pre_E(:, betatime:(betatime + 30*sampling_freq), ch), 2); % Taking first 30 seconds after appearance of first beta burst.
                
                avg_E(:, fo, ch) = mean(pre_E(:, infusetimes(fo)*sampling_freq:(infusetimes(fo) + 30)*sampling_freq, ch), 2); % Taking first 30 seconds after infusion.
                
                plot(avg_E(:, fo, ch), 'r') % Plotting average.
                
                hold on
                
                figure(fo + 10)
                
                subplot(2, 2, ch+2) % Plotting all 30 seconds.
                
                imagesc(1:30*sampling_freq, 1:size(pre_E, 1), pre_E(:, infusetimes(fo)*sampling_freq:(infusetimes(fo) + 30)*sampling_freq, ch))
                
                axis xy
            
            end
            
            hold on
            
        end
        
    end
    
    mean_E(:, pd, :) = mean(avg_E, 2); % Mean over animals.
    
    std_E(:, pd, :) = std(avg_E, 1, 2)/sqrt(size(avg_E, 2)); % Standard error over animals.
    
end


for ch = 1:no_channels
    
        figure(fo + 1)
        
        subplot(1, no_channels, ch)
        
        boundedline(1:length(mean_E), reshape(mean_E(:, ch, :), [no_freqs 2]), std_E(:, ch, :))

        title(chan_labels(ch))
        
        % ylabel(sprintf('log (Power)\n[z-score]'))
        ylabel(sprintf('Power\n(normalized to pre-infusion)'))
        
        xlabel('Frequency (Hz)')
        
        if ch==2
        
            legend('Pre-Infusion','Post-Infusion')
        
        end
        
 end

 % figdir='C:\Users\Administrator\SkyDrive\Boston University\Han Lab\Results PPTs\Carb Paper Figures\';
 
 save_as_pdf(gcf, [subject_mat(1:end-4),'_powerfreq_v2'])