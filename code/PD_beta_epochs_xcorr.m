function PD_beta_epochs_xcorr(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, suffix, ~)

% Default suffix.
if nargin < 6
    
    suffix = '';
    
end

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

xcorr_length = 2*win_size + 1;

data_length = win_size + 1;

norm = data_length - abs((1:xcorr_length) - data_length);

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
            
            All_xcorr = nan(no_epochs, 2*win_size + 1);
            
            parfor e = 1:no_epochs
                
                data_name = beta_list{e}; data_name = [data_name(1:(end - 4)), suffix, '.txt']; 
                
                data = load(data_name);
                
                data_xcorr = xcorr(data(:, 1), data(:, 2), 'coeff');
                
                All_xcorr(e, :) = data_xcorr' ./ norm;
                
            end
            
            save([beta_listname(1:end-5),'_xcorr', suffix, '.mat'], 'All_xcorr')
            
        end
        
    end
    
end