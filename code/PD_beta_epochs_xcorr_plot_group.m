function PD_beta_epochs_xcorr_plot_group(subjects_mat, ~, ~, win_size, ~, suffix, suffix_name)

% Default suffix.
if nargin < 6
    
    suffix = ''; suffix_name = 'LFP';
    
end

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

t = (1:(2*win_size + 1)) - win_size;

for ch = 1:no_channels
    
    figure
    
    channel_coh = nan(2*win_size + 1, 2, 2);
    
    for pd = 1:2
        
        listname = ['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch}, '_', pd_label{pd}, '_xcorr', suffix];
        
        All_xcorr = load([listname, '.txt']);
        
        no_epochs = size(All_xcorr, 1);
        
        channel_coh(:, pd, 1) = nanmean(All_xcorr)';
        
        channel_coh(:, pd, 2) = nanstd(All_xcorr)'/sqrt(no_epochs);
        
    end
        
    mean_data = channel_coh(:, :, 1);
    
    std_data = reshape(channel_coh(:, :, 2), [size(channel_coh, 1) 1 2]);
    std_data = repmat(std_data, [1 2 1]);
    
    boundedline(t', mean_data, std_data)
    
    legend(period_label)
    
    axis tight
    
    title({[chan_labels{ch}, ' High Beta Blocks'];['Mean Cross-Correlation of ', suffix_name]})
    
    save_as_pdf(gcf, listname)
    
end