function MM_beta_epochs_coh_mtm(filenames, chan_labels, infusion_times, outlier_lim, sd_lim, win_size, smooth_size, tbw)

% Computes complex cross-spectrum between channels for each epoch, so
% that coherence and phase of coherence can be calculated. Utilizes
% multitaper method to compute the complex Fourier transform used to
% calculate the cross-spectrum.
%
% SAMPLE CALL:
% MM_beta_epochs_coh_mtm_collect({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},7,2,1000,5000,1*2)
%
% INPUTS:
% 'filenames' is a cell of strings, which are the filenames of files 
% containing data for picking beta segments. The data should contain two
% channels, as columns.
% 'chan_labels' is a cell of strings, which are the labels of channels
% inside each file containing data.
% 'infusion_times' is a vector, with length the same as 'filenames', 
% containing the times (in datapoints) at which each data file switches
% from control to Parkinsonian behavior (or the time of carbachol
% infusion).
% 'outlier_lim' is the number of standard deviations beyond which a spike
% in the LFP is considered an outlier.
% 'sd_lim' is the number of standard deviations defining the cutoff of high
% beta power.
% 'win_size' is the minimum duration (in datapoints, so s*sampling_freq) 
% for which beta must be elevated above the cutoff, to be considered a 
% high beta segment.
% 'smooth_size' is the length of time (in datapoints, so s*sampling_freq)
% over which beta power is smoothed before applying the cutoff.
% 'tbw' is the time bandwidth product.

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

pd_label = {'pre', 'post'};

for fi = 1:length(filenames)
    
    filename = filenames{fi};
    
    for ch = 1:no_channels
        
        for pd = 1:2
            
            beta_listname = [filename,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            beta_list = text_read(beta_listname, '%s%*[^\n]');
            
            no_epochs = length(beta_list);
            
            All_coh = nan(no_epochs, win_size + 1);
            
            parfor e = 1:no_epochs
                
                data_name = beta_list{e};
                
                data = load(data_name);
                
                data_xspec_norm = xspec_mtm(data(:, 1), data(:, 2), tbw);
                
                All_coh(e, :) = data_xspec_norm';
                
            end
            
            save([beta_listname(1:end-5), '_coh_mtm_', num2str(tbw), 'tbw.mat'], 'All_coh')
            
        end
        
    end
    
end

end

function xc_mtm = xspec_mtm(x, y, tbw)

% Takes x and y data (to be cross-correlated) and time-bandwidth product.

length_x = length(x); length_y = length(y);

if length_x ~= length_y
    
    display('Vectors must be the same length.')
    
    return

else
    
    if size(x, 2) ~= 1, x = x'; end
    
    if size(y, 2) ~= 1, y = y'; end
    
    data = [x y];
    
    dps_seq = dpss(length_x, tbw/2);
    
    xspec_norm_est = nan(size(dps_seq));
    
    for s = 1:tbw
        
        fft_est = fft(data.*repmat(dps_seq(:, s), 1, 2));
       
        xspec_est = fft_est(:, 1) .* conj(fft_est(:, 2));
        
        spec_est = sqrt(fft_est .* conj(fft_est));
        
        xspec_norm_est(:, s) = xspec_est./prod(spec_est, 2);
        
    end
    
    xc_mtm = nanmean(xspec_norm_est');

end

end