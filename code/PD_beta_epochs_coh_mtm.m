function PD_beta_epochs_coh_mtm(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)

% Takes typical arguments, plus time-bandwidth product.

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
            
            All_coh = nan(no_epochs, win_size);
            
            for e = 1:no_epochs
                
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