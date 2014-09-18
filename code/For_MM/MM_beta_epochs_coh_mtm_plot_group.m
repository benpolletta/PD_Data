function MM_beta_epochs_coh_mtm_plot_group(filenames, ~, chan_labels, ~, ~, ~, win_size, ~, tbw)

% Plots magnitude (coherence) and phase (phase of coherence) of the expected complex
% cross-spectrum between channels. Utilizes jackknife resampling to compute
% standard deviation of coherence and phase of coherence.
%
% SAMPLE CALL:
% MM_beta_epochs_coh_mtm_plot_group({'file1.txt','file2.txt'},1000,{'Striatum','Motor
% Ctx.'},7,2,1000,5000,1*2)
%
% INPUTS:
% 'filenames' is a cell of strings, which are the filenames of files 
% containing data for picking beta segments. The data should contain two
% channels, as columns.
% 'chan_labels' is a cell of strings, which are the labels of channels
% inside each file containing data.
% 'infusion_times' is a vector, with length the same as 'filenames', 
% containing the times (in datapoints) at which each data file switches
% from control to Parkinsonian behavior (or the time of carbachol
% infusion).
% 'outlier_lim' is the number of standard deviations beyond which a spike
% in the LFP is considered an outlier.
% 'sd_lim' is the number of standard deviations defining the cutoff of high
% beta power.
% 'win_size' is the minimum duration (in datapoints, so s*sampling_freq) 
% for which beta must be elevated above the cutoff, to be considered a 
% high beta segment.
% 'smooth_size' is the length of time (in datapoints, so s*sampling_freq)
% over which beta power is smoothed before applying the cutoff.
% 'tbw' is the time bandwidth product.

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

% Creating name to call the analysis by - either name of single file in
% 'filenames', or 'All'. NB: will be over-written each time you run a
% multi-file analysis, so should be renamed after running.

if length(filenames) == 1
    
    file_label = filenames{1};
    
else

    file_label = 'All';
    
end

% Labeling periods (pre- and post-infusion).

period_label = {'Pre-Injection','Post-Infusion'};
pd_label = {'pre', 'post'};

% Frequency information.

f = 1000*(0:(win_size - 1))/win_size;

f_indices = f <= 32 & f >= 8;

% Labels for coherence measures (magnitude and phase of expected
% cross-correlation).

m_label = {'coh','phase'};
measure_label = {'Coherence','Phase of Coherence'};
no_measures = length(measure_label);

for ch = 1:4 %no_channels
    
    figure
    
    channel_coh = nan(win_size, 2, 2, 2);
    
    for pd = 1:2
        
        listname = [file_label, '_', ch_label{ch}, '_', pd_label{pd}, '_coh_mtm_', num2str(tbw), 'tbw'];
        
        All_coh_r = load([listname, '_r.txt']);
        
        All_coh_i = load([listname, '_i.txt']);
        
        All_coh = All_coh_r + sqrt(-1)*All_coh_i;
        
        no_epochs = size(All_coh, 1);
        
        if no_epochs > 0
            
            if no_epochs > 2
                
                All_jack = jackknife(@nanmean, All_coh);
                
            else
                
                All_jack = jackknife(@nanmean, [All_coh; nan(3 - no_epochs, size(All_coh, 2))]);
                
            end
            
            channel_coh(:, pd, 1, 1) = abs(nanmean(All_coh))';
            
            channel_coh(:, pd, 2, 1) = sqrt((no_epochs - 1)*nanstd(abs(All_jack)).^2); %/sqrt(no_epochs)';
            
            channel_coh(:, pd, 1, 2) = angle(nanmean(All_coh))';
            
            channel_coh(:, pd, 2, 2) = sqrt((no_epochs - 1)*circ_std(angle(All_jack)).^2); %/sqrt(no_epochs)';
            
        end
        
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
        
        save([listname, '_', m_label{measure}], 'mean_data', 'std_data')
        
    end
    
    save_as_pdf(gcf, listname)
    
end