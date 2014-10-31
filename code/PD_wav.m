function PD_wav(subjects_mat)

load(subjects_mat)

sampling_freq = 1000;

freqs = 1:200;
no_freqs = length(freqs);

bands = [1 4; 4 8; 8 30; 30 100; 100 120; 0 200]; no_bands = size(bands, 1);

no_cycles = linspace(3, 21, no_freqs);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', sampling_freq);

segment_length = size(wavelets, 2);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    base_index = basetime*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    % if isempty(dir([subj_name, '_wt.mat']))
        
        load([subj_name,'_all_channel_data_dec.mat'])
        
        no_secs = min(basetime + 1500, length(PD_dec)/sampling_freq);
        
        t = (1:no_secs*sampling_freq)/sampling_freq;
        
        clear Spec Spec_pct BP BP_norm BP_pct
        
        [Spec, Spec_pct] = deal(nan(no_secs*sampling_freq, no_freqs, 2));
        
        [BP, BP_norm, BP_pct] = deal(nan(no_secs*sampling_freq, no_bands, 2));
        
        for ch = 1:2
            
            data = PD_dec(1:(no_secs*sampling_freq), ch);
            
            data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))];
            
            for f = 1:no_freqs
                
                conv_prod = conv(data_reflected, wavelets(f,:), 'same');
                
                Spec(:, f, ch) = conv_prod((segment_length + 1):(end - segment_length));
                
            end
            
            %% Band power.
            
            for b = 1:no_bands
                
                band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2);
                
                BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2));
                
            end
            
            %% Normalize by total power.
            
            BP_norm(:, :, ch) = BP(:, :, ch)./repmat(BP(:, end, ch), 1, no_bands);
            
            Spec_norm(:, :, ch) = Spec(:, :, ch)./repmat(BP(:, end, ch), 1, no_freqs);
            
            %% Baseline normalize.
            
            baseline_mean = mean(abs(Spec_norm(t <= basetime, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
            
            Spec_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
            
            baseline_BP = mean(BP_norm(t <= basetime, :, ch));
            
            BP_pct(:, :, ch) = 100*BP_norm(:, :, ch)./ones(size(BP_norm(:, :, ch)))*diag(baseline_BP) - 100;
            
            % figure(1)
            % 
            % subplot(2, 1, ch)
            % 
            % imagesc(t - basetime, freqs, zscore(abs(Spec(:, :, ch)))')
            % 
            % cl = caxis; caxis([cl(1) .25*cl(2)])
            % 
            % axis xy
            % 
            % xlabel('Time (sec)'); ylabel('Hz');
            % 
            % title(['Gabor Spectrogram of ', folder, ', ', chan_labels{ch}])
            % 
            % grid on
            % 
            % figure(2)
            % 
            % subplot(2, 1, ch)
            % 
            % imagesc(t - basetime, freqs, Spec_pct(:, :, ch)')
            % 
            % cl = caxis; caxis([cl(1) .25*cl(2)])
            % 
            % axis xy
            % 
            % xlabel('Time (sec)'); ylabel('Hz');
            % 
            % title(['Baseline Normalized Gabor Spectrogram of ', folder, ', ', chan_labels{ch}])
            % 
            % grid on
            
        end
        
        save([subj_name, '_wt.mat'], '-v7.3', 'Spec', 'Spec_pct', 't', 'basetime', 'freqs', 'BP', 'BP_norm', 'BP_pct')
        
        % suffix = {'', '_pct'};
        %
        % for s = 1:length(suffix)
        % 
        %     try
        % 
        %         save_as_pdf(s, [subj_name, '_wt', suffix{s}])
        % 
        %     catch
        % 
        %         saveas(s, [subj_name, '_wt', suffix{s}, '.fig'])
        % 
        %     end
        % 
        % end
        
    %end
    
end



% %% Morlet wavelet spectrogram.
%
% MorletFourierFactor = 4*pi/(6+sqrt(2+6^2));
% freq = 1:200;
% scales = 1./(freq*MorletFourierFactor);
%
% PD_sig =  struct('val',PD_dec,'period',1/sampling_freq);
% cwtPD = cwtft(PD_sig,'scales',scales);
% scales = cwtPD.scales;
% freq = 1./(scales*MorletFourierFactor);
%
% figure;imagesc(t,freq,zscore(abs(cwtPD.cfs)));set(gca,'YDir','normal');
% xlabel('time (sec)'); ylabel('Pseudo-frequency');
% title('Morlet spectrogram of PD data')
% set(gca,'YTick',freq(1:10:length(freq)))
% grid on
% ylim([freq(1) freq(end)])

% %% Gabor spectrogram.
% 
% NFFT    = 2^12;
% 
% tw      = -NFFT/2+1:NFFT/2;
% 
% sigma   = .2; %[sec]
% 
% sigSamp = sigma*sampling_freq;
% 
% w       = sqrt(sqrt(2)/sigSamp)*exp(-pi*tw.*tw/(sigSamp^2)); % Gaussian window.
% 
% overlap = NFFT-1;
% 
% [SpecPD, T, F] = spectrogram(PD_dec(1:(base_index + 1500*sampling_freq), ch), w, overlap, NFFT, sampling_freq); % Derive gaussian windowed spectrogram.