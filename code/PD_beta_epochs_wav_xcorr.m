function PD_beta_epochs_wav_xcorr(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, freqs)

if nargin < 6
    
    freqs = 8:2:32;
    
end

no_freqs = length(freqs);

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    for ch = 1:no_channels
        
        for pd = 1:2
            
            beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            beta_list = text_read(beta_listname, '%s%*[^\n]');
            
            no_epochs = length(beta_list);
            
            All_xcorr = nan(no_epochs, (2*win_size + 1)*no_freqs, 2);
            
            for e = 1:no_epochs
                
                data_name = beta_list{e};
                
                data = load(data_name);
                
                data_wav_xcorr = wav_xcorr(data, freqs, 1000);
                
                All_xcorr(e, :, :) = reshape(data_wav_xcorr, 1, (2*win_size + 1)*no_freqs, 2);
                
            end
            
            save([beta_listname(1:end-5),'_wav_xcorr.mat'], 'All_xcorr')
            
        end
        
    end
    
end

end

function wxc = wav_xcorr(data, freqs, sampling_freq)

    no_freqs = length(freqs);
    
    data_length = size(data, 1);
    
    xcorr_length = 2*data_length - 1;
        
    norm = data_length - abs((1:xcorr_length) - data_length);
    
    norm = norm / max(norm);
    
    wxc = nan(xcorr_length, no_freqs, 2);
    
    wav = dftfilt3(freqs, 7, sampling_freq);
    
    for f = 1:no_freqs
       
        x_wav = conv(data(:, 1), wav{f}', 'same');
        
        y_wav = conv(data(:, 2), wav{f}', 'same');
        
        real_xcorr = xcorr(real(x_wav), real(y_wav), 'coeff');
        
        wxc(:, f, 1) = real_xcorr ./ norm'; % ./ (denom_1 .* denom_2);
        
        abs_xcorr = xcorr(abs(x_wav), abs(y_wav), 'coeff');
        
        wxc(:, f, 2) = abs_xcorr ./ norm';
        
    end

end
        
% denom_1 = sqrt(xcorr(x_wav.^2, ones(size(y_wav))));
%
% denom_2 = sqrt(xcorr(ones(size(x_wav)), y_wav.^2));