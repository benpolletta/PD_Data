function PD_beta_epochs_coh_mtm_collect(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, tbw)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

for ch = 1:no_channels
    
    for pd = 1:2
        
        fid1 = fopen(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch},...
            '_', pd_label{pd}, '_', par_name, '_coh_mtm_', num2str(tbw), 'tbw_r.txt'], 'w');
        
        fid2 = fopen(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch},...
            '_', pd_label{pd}, '_', par_name, '_coh_mtm_', num2str(tbw), 'tbw_i.txt'], 'w');
        
        for fo = 1:length(folders)
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            subj_name = [folder,'/',prefix];
            
            beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            load([beta_listname(1:end-5),'_coh_mtm_', num2str(tbw), 'tbw.mat'])
            
            format = make_format(size(All_coh, 2), 'f');
            
            fprintf(fid1, format, real(All_coh)');
            
            fprintf(fid2, format, imag(All_coh)');
            
        end
        
    end
    
end