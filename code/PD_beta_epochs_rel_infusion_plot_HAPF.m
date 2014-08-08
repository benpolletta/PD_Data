function PD_beta_epochs_rel_infusion_plot_HAPF(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

sampling_freq = 1000;

smooth_winsize = 50;

f_bins = 8:4:32;

f_labels = textscan(num2str(f_bins), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

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
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    for ch = 1:2%no_channels
        
        beta_name = [subj_name,'_',ch_label{ch},'_beta'];
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            beta_pbf_name = [beta_listname(1:end-5),'_pbf_dp.txt'];
            
            % if isempty(dir(beta_pbf_name))
                
                fid = fopen(beta_pbf_name, 'w');
                
                no_blocks = max(blocks);
                
                %% Computing smoothed frequency and phase difference for each block.
                
                for b = 1:no_blocks
                    
                    beta_start = min(beta_starts(blocks == b));
                    
                    beta_end = max(beta_ends(blocks == b));
                    
                    t = (beta_start:beta_end)/sampling_freq;
                    
                    LFP_block = PD_dec(beta_start:beta_end, :);
                    
                    H_block = H(beta_start:beta_end, :, 3);
                    
                    A_block = A(beta_start:beta_end, :, 3);
                    
                    P_block = unwrap(P(beta_start:beta_end, :, 3));
                    
                    P_diff = -diff(P_block,[],2);
                    
                    P_diff = angle(exp(sqrt(-1)*P_diff));
                    
                    Pd_flipped = [flipud(P_diff(1:smooth_winsize)); P_diff; flipud(P_diff((end-smooth_winsize+1):end))];
                    
                    Pd_conv = angle(conv(exp(sqrt(-1)*Pd_flipped), hann(smooth_winsize)/sum(hann(smooth_winsize)), 'same'));
                    
                    Pd_smooth = Pd_conv((smooth_winsize + 1):(end - smooth_winsize));
                    
                    F = diff(P_block)/(2*pi*(1/sampling_freq));
                    
                    F_smooth = nan(size(Pd_smooth,1),2);
                    
                    for ch1 = 1:2
                        
                        F_flipped = [flipud(F(1:smooth_winsize,ch1)); F(:,ch1); flipud(F((end-smooth_winsize+1):end,ch1))];
                        
                        % F_conv = conv(F_flipped,hann(smooth_winsize)/sum(hann(smooth_winsize)),'same');
                        F_conv = conv(F_flipped,ones(smooth_winsize,1)/smooth_winsize,'same');
                        
                        F_smooth(:,ch1) = F_conv(smooth_winsize + (1:size(Pd_smooth,1)));
                        
                    end
                    
                    figure;
                    
                    %% Plot LFP.
                    
                    subplot(5, 1, 1)
                    
                    plot(t, LFP_block)
                    
                    title([folder, ', ', period_label{pd}, ', Block ', num2str(b), ', Raw LFP'])
                    
                    box off
                    
                    axis tight
                    
                    legend(chan_labels, 'Location', 'NorthWest')
                    
                    %% Plot Bandpassed Oscillation and Phase Diff.
                    
                    subplot(5, 1, 2)
                    
                    [ax, ~, h2] = plotyy(t, H_block, t, Pd_smooth);
                    
                    box off
                    
                    title(['\beta Oscillation and Phase Difference (',chan_labels{1},'-',chan_labels{2},')'])
                    
                    axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
                    
                    set(ax,{'ycolor'},{'black';'black'})
                    
                    set(get(ax(1),'YLabel'),'String','\beta Oscillation')
                    
                    set(get(ax(2),'YLabel'),'String','\Delta \phi')
                    
                    %% Plot Amplitude and Phase Diff.
                    
                    subplot(5, 1, 3)
                    
                    [ax, ~, h2] = plotyy(t, A_block, t, Pd_smooth);
                    
                    box off
                    
                    title(['\beta Amp. and Phase Difference (',chan_labels{1},'-',chan_labels{2},')'])
                    
                    axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
                    
                    set(ax,{'ycolor'},{'black';'black'})
                    
                    set(get(ax(1),'YLabel'),'String','\beta Amplitude')
                    
                    set(get(ax(2),'YLabel'),'String','\Delta \phi')
                    
                    %% Plot Frequency and Phase-Locking.
                    
                    subplot(5, 1, 4)
                    
                    [ax, ~, h2] = plotyy(t, F_smooth, t, Pd_smooth);
                    
                    box off
                    
                    title(['Instantaneous Freq. and Phase Difference (',chan_labels{1},'-',chan_labels{2},')'])
                    
                    xlabel('Time (s)')
                    
                    axis(ax(1),'tight'), axis(ax(2),'tight'), xlim(ax(2),xlim(ax(1)))
                    
                    set(ax,{'ycolor'},{'black';'black'})
                    
                    set(get(ax(1),'YLabel'),'String','\beta Frequency')
                    
                    set(get(ax(2),'YLabel'),'String','\Delta \phi')
                    
                    
                    %% Plot Roseplots.
                    
                    for ch1 = 1:2
                        
                        subplot(5, 4, 4*4 + ch1)
                        
                        [MRV_vec, ~, ~, conf_vec] = rose_plot(Pd_smooth, F_smooth(:,ch1), 20, f_bins);
                        
                        title(['\Delta \phi by ',chan_labels{ch1},' Freq.'])
                        
                    end
                    
                    
                    freezeColors
                    
                    %% Plot Colorplots.
                    
                    colormap('jet')
                    
                    for ch1 = 1:2
                        
                        subplot(5, 4, 4*4 + 2 + ch1)
                        
                        [histogram, bins] = hist3([Pd_smooth F_smooth(:,ch1)], [50 50]);
                        
                        imagesc(bins{2}, bins{1}, histogram)
                        
                        axis xy
                        
                        xlim([10 30])
                        
                        xlabel('Freq. (Hz)')
                        
                        ylabel('\Delta \phi (rad)')
                        
                        title(['\Delta \phi by ',chan_labels{ch1},' Freq.'])
                        
                    end
                    
                    
                    %%
                    
                    save_as_pdf(gcf,[beta_name,'_',pd_label{pd},'_',par_name,'_block',num2str(b)])
                    
                end
                
                close('all')
                
            % end
            
        end
        
    end
    
end