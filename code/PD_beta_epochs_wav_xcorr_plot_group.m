function PD_beta_epochs_wav_xcorr_plot_group(subjects_mat, ~, ~, win_size, ~, freqs)

color_map = colormap('cool'); color_map = flipud(color_map);

% Default suffix.
if nargin < 6
    
    freqs = 8:2:32;
    
end

no_freqs = length(freqs);

c_order = [linspace(1,0,no_freqs); abs(linspace(1,0,no_freqs)-.5); linspace(0,1,no_freqs)]';

load(subjects_mat)

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

period_label = {'Pre-Injection','Post-Injection'};
pd_label = {'pre', 'post'};

measure_label = {'Beta Osc.','Beta Amp.'};

t = (1:(2*win_size + 1)) - win_size;

for ch = 1:no_channels
    
    figure
    
    for pd = 1:2
        
        clear All_xcorr
        
        listname = ['All_', subjects_mat(1:(end-length('_subjects.mat'))), '_', ch_label{ch}, '_', pd_label{pd}, '_wav_xcorr'];
        
        All_xcorr(:, :, 1) = load([listname, '_h.txt']);
        
        All_xcorr(:, :, 2) = load([listname, '_a.txt']);
        
        no_epochs = size(All_xcorr, 1);
        
        mean_xcorr = reshape(nanmean(All_xcorr), size(All_xcorr, 2)/no_freqs, no_freqs, 2);
    
        std_xcorr = reshape(nanstd(All_xcorr), size(All_xcorr, 2)/no_freqs, no_freqs, 2)/sqrt(no_epochs);
        
        for measure = 1:2
           
            subplot(2, 2, (measure - 1)*2 + pd)
            
            std_data = std_xcorr(:, :, measure);
            
            std_data = reshape(std_data, size(std_data, 1), 1, size(std_data, 2));
            
            std_data = repmat(std_data, [1 2 1]); 
            
            boundedline(t', mean_xcorr(:, :, measure), std_data, 'cmap', c_order)
            
            colormap(c_order)
            
            colorbar('YTick', ((0:no_freqs) + .5)/no_freqs, 'YTickLabel', freqs)
            
            axis tight
            
            title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Mean Cross-Corr., ', measure_label{measure}]})
            
        end
        
    end
    
    save_as_pdf(gcf, listname)
    
end