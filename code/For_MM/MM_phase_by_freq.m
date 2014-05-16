function MM_phase_by_freq(filenames,channel_labels,sampling_freq,no_p_bins,f_bins)

% INPUTS:
% 'filenames' - cell of strings, each the name of a file containing data for which 
% a rose plot should be calculated.
% 'channel_labels' - cell arrary of strings, containing labels for each row
% in each file (e.g., {'Striatum','Motor Ctx.'}).
% 'sampling freq' - integer, sampling frequency of the data.
% 'no_p_bins' - number of phase bins to use for each rose plot.
% 'f_bins' - vector containing edges of frequency bins. A separate rose
% plot will be calculated for all timepoints having instantaneous frequency
% in each frequency bin.
%
% SAMPLE CALL:
% MM_phase_by_freq({'test.mat'},{'Motor Ctx.','Striatum'},200,24,6:4:34)

for f = 1:length(filenames)
    
    filename = filenames{f};
    
    phase_name = [filename,'_beta_PbyF'];
    
    data = load(filename);
    
    if isstruct(data)
        
        fields = fieldnames(data);
        
        data = getfield(data,fields{1});
        
    end
    
    no_channels = size(data,2);
    
    BP = eegfilt(data',sampling_freq,10,30);
    
    H = hilbert(BP'); A = abs(H); P = angle(H);
    
    P_unwrapped = unwrap(P);
    
    P_diff = angle(exp(sqrt(-1)*(diff(P_unwrapped,[],2))));
    
    F = diff(P_unwrapped)/(2*pi*(1/sampling_freq));
    
    F_smooth = nan(size(F)); P_diff_smooth = nan(size(P));
    
    for ch = 1:2
       
        F_smooth(:,ch) = conv(F(:,ch),ones(50,1)/50,'same');
        
    end
       
    P_diff_smooth = angle(conv(exp(sqrt(-1)*P_diff),ones(500,1)/500,'same'))/pi;
    
    figure;
    
    subplot(2,1,1)
    
    plot((1:size(F,1))'/sampling_freq,F_smooth)
    
    axis tight
    
    title(filename)
    
    xlabel('Time (s)')
    
    ylabel('Frequency (Hz)')
    
    legend({'Striatal','Motor'})
    
    subplot(2,1,2)
    
    plot((1:size(P,1))'/sampling_freq,P_diff_smooth,'k')
    
    axis tight
    
    xlabel('Time (s)')
    
    ylabel({'Phase Difference';['(',channel_labels{2},' - ',channel_labels{1},', \pi)']})
    
    MR_vec = nan(2, length(f_bins)-1);
    
    p_hist = nan(no_p_bins+1, length(f_bins)-1, 2);
    
    freqs = nan(2, length(f_bins)-1);
    
    figure;
    
    for ch = 1:2
        
        subplot(1,2,ch)
        
        [MR_vec(ch,:), p_hist(:,:,ch), freqs(ch,:)] = rose_plot(P_diff, F_smooth(:,ch), no_p_bins, f_bins);
        
        title({[filename,', Phase Difference (',channel_labels{2},' - ',channel_labels{1},')'];['By ',channel_labels{ch},' Frequency']})
        
    end
    
    save([phase_name,'.mat'],'H','A','P','P_diff','F','F_smooth','MR_vec','p_hist','freqs')
    
    save_as_pdf(gcf,phase_name)
    
end