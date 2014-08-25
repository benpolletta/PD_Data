function PD_beta_epochs_xcorr_plot_individual(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size, suffix, suffix_name)

% Default suffix.
if nargin < 6
    
    suffix = ''; suffix_name = 'LFP';
    
end

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

t = (1:(2*win_size + 1))' - win_size;

for ch = 1:no_channels
    
    for fo = 1:length(folders)
        
        figure
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        subj_name = [folder,'/',prefix];
        
        subject_coh = nan(2*win_size + 1, 2, 2);
        
        for pd = 1:2
            
            beta_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'.list'];
            % beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
            
            load([beta_listname(1:end-5),'_xcorr',suffix,'.mat'])
            
            no_epochs = size(All_xcorr, 1);
            
            if no_epochs > 0
                
                subplot(1, 3, pd)
                
                plot(t, All_xcorr')
                
                xlim([min(t) max(t)])
                
                if pd == 1
                    
                    ylabel('Cross Correlation')
                    
                end
                
                xlabel('Time Lag (ms)')
                
                title({[folder, ', ', chan_labels{ch}, ' High Beta Blocks'];['Cross-Correlation of ',suffix_name,', ', period_label{pd}]})
                
                if no_epochs > 1
                    
                    subject_coh(:, pd, 1) = nanmean(All_xcorr)';
                    
                    subject_coh(:, pd, 2) = nanstd(All_xcorr)'/sqrt(no_epochs);
                    
                else
                    
                    subject_coh(:, pd, 1) = All_xcorr';
                    
                end
                
            end
            
        end
        
        mean_data = subject_coh(:, :, 1);
        
        std_data = reshape(subject_coh(:, :, 2), [size(subject_coh, 1) 1 2]);
        std_data = repmat(std_data, [1 2 1]);
        
        subplot(1, 3, 3)
        
        boundedline(t, mean_data, std_data)
        
        legend(period_label)
        
        axis tight
        
        xlabel('Time Lag (ms)')
        
        title({[folder, ', ', chan_labels{ch}, ' High Beta Blocks'];'Mean Cross-Correlation'})
        
        save_as_pdf(gcf, [subj_name,'_',ch_label{ch},'_beta_',par_name,'_xcorr',suffix])
        
    end
    
end

end