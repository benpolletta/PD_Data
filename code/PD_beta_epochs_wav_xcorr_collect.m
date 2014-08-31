function PD_beta_epochs_wav_xcorr_collect(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, ~)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

measure_label = {'h', 'a'};

for ch = 1:no_channels
    
    for pd = 1:2
        
        for m = 1:2
            
            fid_vec(m) = fopen(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch},...
                '_', pd_label{pd}, '_wav_xcorr_', measure_label{m}, '.txt'], 'w');
            
        end
        
        for fo = 1:length(folders)
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            subj_name = [folder,'/',prefix];
            
            beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            load([beta_listname(1:end-5), '_wav_xcorr.mat'])
            
            format = make_format(size(All_xcorr, 2), 'f');
            
            for m = 1:2
            
                fprintf(fid_vec(m), format, All_xcorr(:, :, m)');
            
            end
                
        end
        
        fclose('all');
        
    end
    
end