function PD_beta_epochs_coh_plot_group(subjects_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

f = 1000*(0:win_size)/win_size;

f_indices = f <= 32 & f >= 8;

measure_label = {'Coherence','Phase of Coherence'};
no_measures = length(measure_label);

for ch = 1:no_channels
    
    figure
    
    channel_coh = nan(win_size + 1, 2, 2, 2);
    
    for pd = 1:2
        
        listname = ['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch}, '_', pd_label{pd}, '_coh'];
        
        All_coh_r = load([listname, '_r.txt']);
        
        All_coh_i = load([listname, '_i.txt']);
        
        All_coh = All_coh_r + sqrt(-1)*All_coh_i;
        
        All_jack = jackknife(@nanmean, All_coh);
        
        no_epochs = size(All_coh, 1);
        
        channel_coh(:, pd, 1, 1) = abs(nanmean(All_coh))';
        
        channel_coh(:, pd, 2, 1) = sqrt((no_epochs - 1)*nanstd(abs(All_jack)).^2); %/sqrt(no_epochs)';
        
        channel_coh(:, pd, 1, 2) = angle(nanmean(All_coh))';
        
        channel_coh(:, pd, 2, 2) = sqrt((no_epochs - 1)*circ_std(angle(All_jack)).^2); %/sqrt(no_epochs)';
        
    end
    
    for measure = 1:no_measures
        
        subplot(1, 2, measure)
        
        mean_data = channel_coh(:, :, 1, measure);
        
        std_data = reshape(channel_coh(:, :, 2, measure), [size(channel_coh, 1) 1 2]);
        std_data = repmat(std_data, [1 2 1]);
        
        boundedline(f(f_indices)', mean_data(f_indices, :), std_data(f_indices, :, :))
        
        legend(period_label)
        
        axis tight
        
        title({[chan_labels{ch}, ' High Beta Blocks'];[measure_label{measure}]})
        
    end
    
    save_as_pdf(gcf, listname)
    
end