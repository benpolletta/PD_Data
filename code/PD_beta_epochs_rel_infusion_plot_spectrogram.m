function PD_beta_epochs_rel_infusion_plot_spectrogram(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

sampling_freq = 1000;

smooth_winsize = 50;

freqs = 0:50;
no_freqs = length(freqs);

no_cycles = 7*ones(no_freqs,1);

wavelets = dftfilt3(freqs, no_cycles, sampling_freq, 'winsize', sampling_freq);

segment_length = size(wavelets, 2);

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

% ch_index = {1, 2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

fid_mat = nan(no_channels,1);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
                    
    A_pct = A(:, :, 3)./sum(A, 3);
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    for ch = 1:no_channels
        
        beta_name = [subj_name,'_',ch_label{ch},'_beta'];
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            %beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
            
            % if isempty(dir(beta_pbf_name))
                
                %fid = fopen(beta_pbf_name, 'w');
                
                no_blocks = max(blocks);
                
                %% Computing smoothed frequency and phase difference for each block.
                
                for b = 1:no_blocks
                    
                    %% Retrieving LFP.
                    
                    beta_start = min(beta_starts(blocks == b));
                    
                    beta_end = max(beta_ends(blocks == b));
                    
                    beta_length = beta_end - beta_start + 1;
                    
                    pad_length = max(333, round(.2*beta_length));
                    
                    plot_start = max(beta_start - pad_length, 1);
                    
                    plot_end = min(beta_end + pad_length, length(PD_dec));
                    
                    t = (plot_start:plot_end)/sampling_freq;
                    
                    LFP_block = PD_dec(plot_start:plot_end, :);
                    
                    H_block = H(plot_start:plot_end, :, 3);
                    
                    P_block = P(plot_start:plot_end, :, 3);
                    
                    A_pct_block = A_pct(plot_start:plot_end, :);
                    
                    %% Calculating wavelet spectrogram.
                    
                    wt_temp = zeros(no_freqs, size(LFP_block, 1), 2);
                    
                    for ch1 = 1:2
                        
                        segment = LFP_block(:, ch1);
                        
                        data_length = size(segment, 1);
                        
                        flip_length = min(segment_length, data_length);
                        
                        segment_reflected = [flipud(segment(1:flip_length)); segment; flipud(segment((end - flip_length + 1):end))];
                        
                        for f = 1:no_freqs
                            
                            conv_prod = conv(segment_reflected, wavelets(f,:), 'same');
                            
                            wt_temp(f, :, ch1) = conv_prod((flip_length + 1):(end - flip_length))';
                            
                        end
                        
                    end
                    
                    %% Calculating phase difference.
                    
                    P_diff = -diff(P_block,[],2);
                    
                    P_diff = angle(exp(sqrt(-1)*P_diff));
                    
                    Pd_flipped = [flipud(P_diff(1:smooth_winsize)); P_diff; flipud(P_diff((end-smooth_winsize+1):end))];
                    
                    Pd_conv = angle(conv(exp(sqrt(-1)*Pd_flipped), hann(smooth_winsize)/sum(hann(smooth_winsize)), 'same'));
                    
                    Pd_smooth = Pd_conv((smooth_winsize + 1):(end - smooth_winsize));
                    
                    figure;
                    
                    %% Plot Spectrogram.
                    
                    for ch1 = 1:2
                        
                        subplot(4, 2, 2*ch1 - 1)%(6, 1, 2*ch1 - 1)
                        
                        imagesc(t, freqs, abs(wt_temp(:, :, ch1)))
                        
                        title({['Wavelet Spectrogram, ', chan_labels{ch1}];[chan_labels{ch}, ' High Beta, ', folder, ', ', period_label{pd}, ', Block ', num2str(b)]}) 
                       
                        axis xy
                        
                        hold on, plot(t', ones(length(t), 2)*diag([10 30]), ':w'), plot(repmat([beta_start, beta_end]/sampling_freq, 2, 1), diag([freqs(1) freqs(end)])*ones(2), 'w')
                        
                        subplot(4, 2, 2*ch1)%(6, 1, 2*ch1)
                        
                        imagesc(t, freqs, zscore(abs(wt_temp(:, :, ch1))')')
                        
                        title('Frequency Normalized Wavelet Spectrogram') 
                       
                        axis xy
                        
                        hold on, plot(t', ones(length(t), 2)*diag([10 30]), '--w'), plot(repmat([beta_start, beta_end]/sampling_freq, 2, 1), diag([freqs(1) freqs(end)])*ones(2), 'w')
                        
                    end
                        
                    %% Plot LFP.
                    
                    subplot(4, 1, 3)%(6, 1, 5)
                    
                    [ax, h1, h2] = plotyy(t, LFP_block, t, A_pct_block);
                    
                    hold on, plot(repmat([beta_start, beta_end]/sampling_freq, 2, 1), diag([max(max(LFP_block)) min(min(LFP_block))])*ones(2), 'k')
                    
                    title(['Raw LFP, ', folder, ', ', period_label{pd}, ', Block ', num2str(b)])
                    
                    box off
                    
                    axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1))) %axis tight
                    
                    set(ax,{'ycolor'},{'black';'black'})
                    
                    set(get(ax(1),'YLabel'),'String','LFP')
                    
                    set(get(ax(2),'YLabel'),'String','Amplitude (% \beta)')
                    
                    legend(h1, chan_labels, 'Location', 'NorthWest')
                    
                    legend(h2, chan_labels, 'Location', 'SouthWest')
                    
                    %% Plot Bandpassed Oscillation and Phase Diff.
                    
                    subplot(4, 1, 4)%(6, 1, 6)
                    
                    [ax, ~, ~] = plotyy(t, real(H_block), t, Pd_smooth);
                    
                    hold on, plot(t, zeros(length(t), 1), ':k'), plot(repmat([beta_start, beta_end]/sampling_freq, 2, 1), diag([max(max(real(H_block))) min(min(real(H_block)))])*ones(2), 'k')
                    
                    box off
                    
                    title(['\beta Oscillation and Phase Difference (',chan_labels{1},'-',chan_labels{2},')'])
                    
                    axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
                    
                    set(ax,{'ycolor'},{'black';'black'})
                    
                    set(get(ax(1),'YLabel'),'String','\beta Oscillation')
                    
                    set(get(ax(2),'YLabel'),'String','\Delta \phi')
                    
                    try
                    
                        save_as_pdf(gcf,[beta_name,'_',pd_label{pd},'_',par_name,'_block',num2str(b),'_wt'])
                    
                    catch
                       
                        display([beta_name,'_',pd_label{pd},'_',par_name,'_block',num2str(b),'_wt.pdf could not be saved.'])
                        
                        saveas(gcf,[beta_name,'_',pd_label{pd},'_',par_name,'_block',num2str(b),'_wt.fig'])
                        
                    end
                        
                end
                
                close('all')
                
            % end
            
        end
        
    end
    
end