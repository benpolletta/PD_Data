function PD_beta_epochs_length_histogram(subject_mat, outlier_lim, sd_lim, win_size, smooth_size, no_bins)

if isscalar(win_size)

    par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

elseif length(win_size) == 2
   
    par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size(1)), 'to', num2str(win_size(2)),'win_',num2str(smooth_size),'smooth'];
    
else
    
    display('win_size must be a scalar (lower limit of beta epoch length) or an interval (lower and upper limits).')
    
    return
    
end
    
load(subject_mat)

pd_label = {'pre', 'post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_index = {1, 2, 1:2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);
    
figure;

All_beta_lengths = cell(2, no_channels);

All_beta_length_histograms = nan(no_bins, 2, no_channels);

All_beta_length_bins = nan(no_bins, 2, no_channels);

[r, c] = subplot_size(no_channels);

for ch = 1:no_channels
    
    for pd = 1:2
        
        for fo = 1:length(folders)
            
            folder = folders{fo};
            
            prefix = prefixes{fo};
            
            subj_name = [folder,'/',prefix];
            
            beta_win_listname = [subj_name,'_',ch_label{ch},'_beta_',pd_label{pd},'_',par_name,'_win.list'];
            
            [block, beta_starts, beta_ends] = text_read(beta_win_listname, '%f%f%f%*[^\n]');
            
            if ~isempty(block)
                
                if isscalar(win_size)
                    
                    no_blocks = max(block);
                    
                    for b = 1:no_blocks
                        
                        block_starts(b) = min(beta_starts(block == b));
                        
                        block_ends(b) = max(beta_ends(block == b));
                        
                    end
                    
                    beta_lengths = (block_ends - block_starts + 1)';
                    
                else
                    
                    beta_lengths = beta_ends - beta_starts + 1;
                    
                end
                
                All_beta_lengths{pd, ch} = [All_beta_lengths{pd, ch}; beta_lengths];
                
            end
            
        end
        
        max_beta_length = max(max(All_beta_lengths{pd, ch}));
        
        min_beta_length = min(min(All_beta_lengths{pd, ch}));
        
        if max_beta_length > 0
            
            bin_edges = linspace(min_beta_length, max_beta_length + 1, no_bins + 1);
            
            bin_centers = (bin_edges(1 : end - 1) + bin_edges(2 : end)) / 2;
            
            [h, b] = histc(All_beta_lengths{pd, ch}, bin_edges);
            
            All_beta_length_histograms(:, pd, ch) = h(1 : end - 1);
            
            All_beta_length_bins(:, pd, ch) = bin_centers;
            
        end
        
    end
    
    subplot(r, c, ch)
    
    norm_hist = All_beta_length_histograms(:, :, ch).*(ones(size(All_beta_length_histograms(:, :, ch)))*diag(1./nansum(All_beta_length_histograms(:, :, ch))));
    
    loglog(All_beta_length_bins(:, :, ch), All_beta_length_histograms(:, :, ch)) %norm_hist)
    
    legend(period_label)
    
    xlabel('Epoch Length')
    
    ylabel('Proportion Observed') %Times Observed')
    
    title({[chan_labels{ch}, ' Beta > ', num2str(sd_lim), ' S.D.'];'Hist. of \beta Epoch Length'})
          
end

save([subject_mat(1:(end - length('_subject.mat'))), 'beta_', par_name, '_length_hist.mat'], 'All_beta_lengths', 'All_beta_length_histograms', 'All_beta_length_bins')

save_as_pdf(gcf, [subject_mat(1:(end - length('_subject.mat'))), 'beta_', par_name, '_', num2str(no_bins), 'bins_length_hist'])