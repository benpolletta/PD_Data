function make_spectrogram_figures(group_prefix, subject_no)

load([group_prefix, '_subjects.mat'])

folder = folders{subject_no}; 
prefix = prefixes{subject_no}; 
basetime = basetimes(subject_no); 
striatum_channel = find(strcmp(chan_labels, 'Striatum'));
    
load([group_prefix, '_pct_BP_high_2.5_min_secs_by_STR.mat'])
    
load([folder, '/', prefix, '_all_channel_data_dec.mat'])
    
data = load([folder, '/', prefix, '_wt.mat']);
Spec = data.Spec;

load([folder, '/', prefix, '_peaks.mat'])

s_rate = 500;
    
dec_factor = 100; s_rate_dec = s_rate/dec_factor;

t = ((1:length(Spec))/s_rate - basetime)/60;

starts = reshape(All_bp_max_start(subject_no, :, 3, :), 2, 2);
starts = starts(striatum_channel, :);
starts = (starts/s_rate - basetime)/60;

ends = reshape(All_bp_max_end(subject_no, :, 3, :), 2, 2);
ends = ends(striatum_channel, :);
ends = (ends/s_rate - basetime)/60;

for ch = 1:2
    
    %% Plotting spectrogram for whole data recording.
    
    figure(ch)
    
    Ch_spec = abs(Spec(:, :, ch));
    
    dims = size(Ch_spec);
    
    Spec_dec = nanmean(reshape(Ch_spec', dims(2), dec_factor, dims(1)/dec_factor), 2);
    Spec_dec = permute(Spec_dec, [1 3 2]);
    
    t_dec = ((1:length(Spec_dec))/s_rate_dec - basetime)/60;
        
    if strcmp(folder, '130328')
        
        Spec_dec = Spec_dec(1:80, :);
        
    else
        
        Spec_dec = Spec_dec(1:80, t_dec <= 30);
        t_dec(t_dec > 30) = [];
        
    end
    
    subplot(3, 1, 1) % figure, subplot(2, 1, 1) % subplot(4, 1, 1)
    
    imagesc(t_dec, 1:80, Spec_dec(1:80, :))
    
    title([folder, ', Channel ', num2str(ch)])
    
    ylabel('Frequency (Hz)')
    
    % caxis([0 10])
    
    axis xy
    
    hold on
    
    plot([0 0], [0 80], ':w', 'LineWidth', 2)
    
    plot([t_dec(1) t_dec(end)], [30 30], ':w', 'LineWidth', 2)
    
    plot([t_dec(1) t_dec(end)], [8 8], ':w', 'LineWidth', 2)
    
    colors = {'g', 'r'};
    
    for pd = 1:2
        
        plot([1 1]*starts(pd), [0 80], colors{pd})
        plot([1 1]*ends(pd), [0 80], colors{pd})
        
    end
    
%     %% Plotting LFP for whole data recording.
%     
%     subplot(2, 1, 2) % subplot(4, 1, 2)
%     
%     plot(t(t <= 30), PD_dec(t <= 30, ch), 'k')
%     
%     box off
%     
%     axis tight
%     
%     xlabel('Time (min.) Rel. Infusion')
%     
%     ylabel('LFP (mV)')
%     
%     try
%         
%         save_as_eps(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch)])
%         
%     catch
%         
%         save(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch), '.fig'])
%         
%     end
    
end

for ch = 1:2
    
    %% Plotting spectrogram & LFP for 10 seconds.
    
    figure(ch)
    
    Ch_spec = abs(Spec(:, :, ch));
    
    pd_labels = {'Pre-Infusion', 'Post-Infusion'};
    
    for pd = 1:2
        
        midpoint = (starts(pd) + ends(pd))/2;
        
        start_time = midpoint - 5/60; end_time = midpoint + 5/60;
        
        t_index = t >= start_time & t <= end_time;
        
        t_interval = t(t_index)*60;
        
        subplot(3, 2, 2 + pd) % figure, subplot(2, 1, 1) % subplot(4, 2, 4 + pd)
        
        imagesc(t_interval, 1:80, Ch_spec(t_index, 1:80)')
        
        title([folder, ', ', pd_labels{pd}])
        
        ylabel('Frequency (Hz)')
        
        % caxis([0 10])
        
        axis xy
        
        hold on
        
        if sum(Spike_indicator(t_index, ch)) > 0
            
            [spike_wav_nans, ~] = indicator_to_nans(Spike_indicator(t_index, ch), 500, 1:200, linspace(3, 21, 200), [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]);
            
            [patch_X, patch_Y] = nans_to_patch(t_interval, 1:90, spike_wav_nans(:, 1:90));
            
            for p = 1:length(patch_X)
                
                p_handle = patch(patch_X{p}, patch_Y{p}, [.8 .8 .8], 'FaceColor', [.8 .8 .8], 'EdgeColor', 'none');
                
                alpha(p_handle, .625)
                
            end
            
        end
        
        plot([t_interval(1) t_interval(end)], [30 30], ':w', 'LineWidth', 2)
        
        plot([t_interval(1) t_interval(end)], [8 8], ':w', 'LineWidth', 2)
        
        subplot(6, 2, 8 + pd) % subplot(4, 1, 3) % subplot(4, 2, 6 + pd)
        
        plot(t_interval, PD_dec(t_index, ch), 'k')
        
        hold on
        
        spikes_w_nans = Spike_indicator(t_index, ch);
        
        spikes_w_nans(spikes_w_nans == 0) = nan;
        
        plot(t_interval, PD_dec(t_index, ch).*spikes_w_nans, 'vk')
        
        box off
        
        axis tight
        
        xlabel('Time (sec.) Rel. Infusion')
        
        ylabel('LFP (mV)')
        
%         try
%             
%             save_as_eps(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch), '_', pd_labels{pd}])
%             
%         catch
%             
%             save(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch), '_', pd_labels{pd} '.fig'])
%             
%         end
        
    end
    
    try
        
        save_as_pdf(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch), 'altogether'])
        
    catch
        
        save(gcf, [folder, '/', folder, '_spec_for_paper_ch', num2str(ch), 'altogether.fig'])
        
    end
    
end
