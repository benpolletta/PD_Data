function PD_beta_epochs_GC_collect(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

pd_label = {'pre', 'post'};

list_suffix = {'', '_A'};

no_lists = length(list_suffix);

dir_labels = {'ch1toch2', 'ch2toch1'};
no_directions = length(dir_labels);

for ch = 1:no_channels
    
    for list = 1:1%no_lists
        
        for pd = 1:2
            
            fid = fopen(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch},...
                '_', pd_label{pd}, '_', par_name, list_suffix{list}, '_GC.txt'], 'w');
            
            for dir = 1:no_directions
                
                fid_spec(dir) = fopen(['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_',...
                    ch_label{ch}, '_', pd_label{pd}, '_', par_name, list_suffix{list}, '_GCspec_', dir_labels{dir}, '.txt'], 'w');
                
            end
            
            for fo = 1:length(folders)
                
                folder = folders{fo};
                
                prefix = prefixes{fo};
                
                subj_name = [folder,'/',prefix];
                
                beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,list_suffix{list},'.list'];
                % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
                
                load([beta_listname(1:end-5),'_GC.mat'])
                
                format = make_format(size(All_GC, 2), 'f');
                
                fprintf(fid, format, All_GC');
                
                spec_format = make_format(size(All_GC_spec, 2), 'f');
                
                for dir = 1:no_directions
                    
                    fprintf(fid_spec(dir), spec_format, All_GC_spec(:, :, dir)');
                    
                end
                
            end
            
        end
        
    end
    
end