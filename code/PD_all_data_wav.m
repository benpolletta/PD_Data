function PD_all_data_wav(subjects_mat, peak_suffix, freqs, no_cycles, bands)

% Computes wavelet spectrogram for decimated data recordings contained in
% the folders which are the entries in the 'folders' variable (a cell of
% strings) in the .mat file names subjects_mat. Saves output in multiple
% .mat files.
%
% INPUTS:
% subjects_mat - name of .mat file containing list of folders.
% freqs - a vector of center frequencies for the wavelet spectrogram;
%  default is 1:200.
% no_cycles - a vector of numbers of cycles each wavelet contains
%  (effectively gives the bandwidth of the spectrogram at each frequency);
%  default is linspace(3, 7, 200).
% bands - a matrix of band limits, over which spectrogram power will be
%  summed; matrix is n by 2, with n the number of bands, the first
%  column containing the start frequency of each band, and the second column
%  containing the end frequency of each band; default is [1 4;4 8;8 30;30
%  100;120 180;0 200].

subjects_struct = load(subjects_mat);

folders = subjects_struct.folders;

prefixes = subjects_struct.prefixes;

basetimes = subjects_struct.basetimes;

parfor fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    wav_inner(folder, prefix, basetime, freqs, no_cycles, bands, peak_suffix) % See below.
    
end

end

% function wav_inner(folder, prefix, basetime, freqs, no_cycles, bands)
% 
% subj_name = [folder,'/',prefix];
% 
% if isempty(freqs) && isempty(no_cycles) && isempty(bands) % Defaults.
%     
%     freqs = 1:200;
%     
%     no_cycles = linspace(3, 21, 200);
%     
%     bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
%     
%     save_name = subj_name;
%     
% else % Name to use when saving non-default output.
% 
%     save_name = sprintf('%s_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', subj_name, freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
% 
% end
%     
% load([subj_name,'_all_channel_data_dec.mat']) % Loading decimated data.
% 
% no_freqs = length(freqs); 
% 
% no_bands = size(bands, 1);
% 
% wavelet_lengths = no_cycles.*(sampling_freq./freqs); % Lengths of wavelets.
% 
% wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', max(sampling_freq, max(wavelet_lengths))); % Wavelets.
% 
% segment_length = sampling_freq;
% 
% t = (1:length(PD_dec))/sampling_freq - basetime;
% 
% clear Spec Spec_norm Spec_pct Spec_norm_pct BP BP_norm BP_pct BP_norm_pct % Clearing relevant variables.
% 
% %% Wavelet spectrogram.
% 
% [Spec, Spec_norm, Spec_pct, Spec_norm_pct] = deal(nan(length(PD_dec), no_freqs, 2)); % Pre-allocating spectrograms, normalized various ways.
% 
% for ch = 1:2
%     
%     data = PD_dec(:, ch);
%     
%     data_reflected = [flipud(data(1:segment_length)); data; flipud(data((end - segment_length + 1):end))]; % Reflecting to avoid edge artifacts.
%     
%     for f = 1:no_freqs
%         
%         conv_prod = conv(data_reflected, wavelets(f,:), 'same'); % Convolving data with wavelet.
%         
%         Spec(:, f, ch) = conv_prod((segment_length + 1):(end - segment_length)); % Assigning to spectrogram.
%         
%     end
%     
%     %% Baseline normalize (percent change).
%     
%     baseline_mean = mean(abs(Spec(t <= 0, :, ch)));
%     
%     Spec_pct(:, :, ch) = 100*abs(Spec(:, :, ch))./(ones(size(Spec(:, :, ch)))*diag(baseline_mean)) - 100;
%     
%     %% Normalize by total power.
%     
%     Spec_norm(:, :, ch) = abs(Spec(:, :, ch))./repmat(sqrt(sum(abs(Spec(:, :, ch)).^2, 2)), 1, no_freqs);
%     
%     %% Baseline normalize percent of total power (percent change).
%     
%     baseline_mean = mean(abs(Spec_norm(t <= 0, :, ch))); %baseline_std = std(abs(Spec(t <= basetime, :, ch)));
%     
%     Spec_norm_pct(:, :, ch) = 100*abs(Spec_norm(:, :, ch))./(ones(size(Spec_norm(:, :, ch)))*diag(baseline_mean)) - 100;
%     
% end
% 
% % Saving separately for speed.
% save([save_name, '_wt.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec', '-v7.3')
% save([save_name, '_wt_norm.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm', '-v7.3')
% save([save_name, '_wt_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_pct', '-v7.3')
% save([save_name, '_wt_norm_pct.mat'], 'sampling_freq', 't', 'basetime', 'freqs', 'Spec_norm_pct', '-v7.3')
%     
% %% Band power.
% 
% [BP, BP_norm, BP_pct, BP_norm_pct] = deal(nan(length(PD_dec), no_bands, 2));
% 
% for ch = 1:2
% 
%     for b = 1:no_bands
%         
%         band_indices = freqs >= bands(b, 1) & freqs <= bands(b, 2); % Frequency indices for band.
%         
%         BP(:, b, ch) = sqrt(sum(abs(Spec(:, band_indices, ch)).^2, 2)); % Sum of squared absolute value.
%         
%         BP_pct(:, b, ch) = sum(Spec_pct(:, band_indices, ch), 2); % Sum.
%         
%         BP_norm(:, b, ch) = sum(Spec_norm(:, band_indices, ch), 2); % Sum.
%         
%         BP_norm_pct(:, b, ch) = sum(Spec_norm_pct(:, band_indices, ch), 2); % Sum.
%         
%     end
%     
% end
% 
% save([save_name, '_wt_BP.mat'], 'sampling_freq', 't', 'basetime', 'bands', 'BP', 'BP_norm', 'BP_pct', 'BP_norm_pct', '-v7.3')
% 
% end