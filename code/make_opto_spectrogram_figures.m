function make_opto_spectrogram_figures(group_prefix, peak_suffix, subject_no, trial_to_zoom, color_lims, zoom_loc)

close('all')

figure

if size(color_lims, 1) == 1
    
    color_lims = repmat(color_lims, 2, 1);
    
end

if isempty(trial_to_zoom), trial_to_zoom = 5; end

load([group_prefix, '_subjects.mat'])

folder = folders{subject_no}; 
prefix = prefixes{subject_no};
subj_name = [folder, '/', prefix];
basetime = basetimes(subject_no); 
striatum_channel = find(strcmp(chan_labels, 'Striatum'));
    
load([folder, '/', prefix, '_all_channel_data_dec.mat'])

data = load([subj_name, '_wt.mat']);
Spec = data.Spec;

load([subj_name, '_wav_laser_artifacts.mat'])

no_laser_ons = ceil(cumsum(laser_transitions)/2);

Spike_indicator = peak_loader([folder, '/', prefix], peak_suffix, size(Spec, 1));

for ch = 1:2,
   
    Spike_laser_indicator(:, ch) = double(logical(Spike_indicator(:, ch)) | logical(laser_transitions));
    
end

s_rate = 500;

dec_factor = 100; s_rate_dec = s_rate/dec_factor;

t = ((1:length(Spec))/s_rate)/60;

first_ten_endtime = t(find(no_laser_ons > 10, 1));

trial_borders = t(find(no_laser_ons > trial_to_zoom - 1, 1)) + [-5 10]/60;

prelaser_transitions = circshift(laser_transitions, -5*sampling_freq);

laseron_times = t(logical(laser_transitions & prelaser_transitions));

laseroff_times = t(xor(laser_transitions, prelaser_transitions) & laser_transitions);

prelaseron_times = t(xor(laser_transitions, prelaser_transitions) & prelaser_transitions);

times = {prelaseron_times, laseron_times, laseroff_times};

time_colors = {'m', 'w', 'c'}; time_LFP_colors = {'m', 'k', 'c'};

for c = 1:3

    times{c}(times{c} > first_ten_endtime) = [];
    
end

figure()

%% Plotting trial-averaged spectrogram.

trial_datapoints = int16((times{1}(2) - times{1}(1))*60*sampling_freq);

smooth_x_datapoints = sampling_freq/25; smooth_y_datapoints = sampling_freq/25;

mg_x = -3:(6/(smooth_x_datapoints - 1)):3; mg_y = -3:(6/(smooth_y_datapoints - 1)):3;

[smooth_x, smooth_y] = meshgrid(mg_x, mg_y);

smooth_matrix = mvnpdf([smooth_x(:) smooth_y(:)], [0 0], eye(2));

smooth_matrix = reshape(smooth_matrix/sum(smooth_matrix), smooth_x_datapoints, smooth_y_datapoints);

% single_spike_indicator = zeros(sampling_freq, 1);
% 
% single_spike_indicator(floor(sampling_freq/2)) = 1;
%         
% [single_spike_wav_nans, ~] = indicator_to_nans(single_spike_indicator, 500, 1:200, linspace(3, 21, 200), [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]);

for ch = 1:2
    
    chan_index = abs(3*mod(3 - striatal_id(subject_no), 2) - ch);
    
    Ch_spec = abs(Spec(:, :, chan_index));
    
    trials_start = find(t == times{1}(1));
    trials_end = find(t == times{1}(end)) - 1;
    
    Trials_spec = Ch_spec(trials_start:trials_end, :);
    
    Trials_spec = reshape(Trials_spec, size(Trials_spec, 1)/trial_datapoints, trial_datapoints, size(Trials_spec, 2));
    
    Mean_spec = reshape(nanmean(Trials_spec), trial_datapoints, size(Trials_spec, 3))';
    
    t_trial = (1:trial_datapoints)/sampling_freq - 5;
    
    % for f = 1:80
    % 
    %     Mean_spec(f, :) = conv(Mean_spec(f, :), single_spike_wav_nans(f, :), 'same');
    % 
    % end
    
    Mean_spec = conv2(Mean_spec(1:80, t_trial <= 10), smooth_matrix, 'same');
    
    % Mean_spec = conv2(Mean_spec(1:80, t_trial <= 10), ones(sampling_freq/25)/((sampling_freq/25)^2), 'same');
    
    subplot(3, 2, ch) % figure, subplot(2, 1, 1) % subplot(4, 1, 1)
    
    imagesc(t_trial(t_trial <= 10), 1:80, Mean_spec)
    
    set(gca, 'FontSize', 16)
    
    caxis(color_lims(chan_index, :))
    
    title(chan_labels{ch}, 'FontSize', 20) % [folder, ', Channel ', num2str(ch)])
    
    ylabel('Frequency (Hz)', 'FontSize', 16)
    
    xlabel('Time (sec.) Rel. Laser On', 'FontSize', 16)
    
    % caxis([0 10])
    
    axis xy
    
    hold on
    
    plot([t_trial(1) t_trial(end)], [30 30], ':w', 'LineWidth', 2)
    
    plot([t_trial(1) t_trial(end)], [8 8], ':w', 'LineWidth', 2)
    
    plot([-5 -5], [0 80], time_colors{1}, 'LineWidth', 2)
    
    plot([0 0], [0 80], time_colors{2}, 'LineWidth', 2)
    
    plot([5 5], [0 80], time_colors{3}, 'LineWidth', 2)

end

% %% Plotting spectrogram for first 10 trials.
% 
% for ch = 1:2
%     
%     chan_index = abs(3*mod(3 - striatal_id(subject_no), 2) - ch);
%     
%     Ch_spec = abs(Spec(:, :, chan_index));
%     
%     dims = size(Ch_spec);
%     
%     Spec_dec = nanmean(reshape(Ch_spec', dims(2), dec_factor, dims(1)/dec_factor), 2);
%     Spec_dec = permute(Spec_dec, [1 3 2]);
%     
%     t_dec = ((1:length(Spec_dec))/s_rate_dec)/60;
%     
%     subplot(3, 2, ch) % figure, subplot(2, 1, 1) % subplot(4, 1, 1)
%     
%     imagesc(t_dec(t_dec < first_ten_endtime), 1:80, Spec_dec(1:80, t_dec < first_ten_endtime))
%     
%     set(gca, 'FontSize', 16)
%     
%     caxis(color_lims(chan_index, :))
%     
%     title(chan_labels{ch}, 'FontSize', 20) % [folder, ', Channel ', num2str(ch)])
%     
%     ylabel('Frequency (Hz)', 'FontSize', 16)
%     
%     xlabel('Time (min.)', 'FontSize', 16)
%     
%     % caxis([0 10])
%     
%     axis xy
%     
%     hold on
%     
%     plot([t_dec(1) t_dec(end)], [30 30], ':w', 'LineWidth', 2)
%     
%     plot([t_dec(1) t_dec(end)], [8 8], ':w', 'LineWidth', 2)
%     
%     for c = 1:3
%         
%         plot(repmat(times{c}, 2, 1), repmat([0; 80], 1, size(times{c}, 2)), time_colors{c}, 'LineWidth', 1.5)
%         
%     end
%     
% end

%% Plotting spectrogram & LFP for 10 seconds (during trial 5).

for ch = 1:2
    
    chan_index = abs(3*mod(3 - striatal_id(subject_no), 2) - ch);
    
    Ch_spec = abs(Spec(:, :, chan_index));
    
    pd_labels = {'Pre-Infusion', 'Post-Infusion'};
    
    start_time = trial_borders(1); end_time = trial_borders(2);
    
    t_index = t >= start_time & t <= end_time;
    
    for c = 1:3
    
        trial_times{c} = times{c}(times{c} >= start_time & times{c} <= end_time)*60;
        
    end
    
    t_interval = t(t_index)*60;
    
    subplot(3, 2, 2 + ch) % figure, subplot(2, 1, 1) % subplot(4, 2, 4 + pd)
    
    imagesc(t_interval, 1:80, Ch_spec(t_index, 1:80)')
    
    set(gca, 'FontSize', 16)
    
    % title(pd_labels{pd}, 'FontSize', 16)
    
    ylabel('Frequency (Hz)', 'FontSize', 16)
    
    caxis(color_lims(chan_index, :)) % caxis([0 10])
    
    axis xy
    
    hold on
    
    if sum(Spike_laser_indicator(t_index, chan_index)) > 0
        
        [spike_wav_nans, ~] = indicator_to_nans(Spike_laser_indicator(t_index, chan_index), 500, 1:200, linspace(3, 21, 200), [1 4; 4 8; 8 30; 30 100; 120 180; 0 200]);
        
        [patch_X, patch_Y] = nans_to_patch(t_interval, 1:90, spike_wav_nans(:, 1:90));
        
        for p = 1:length(patch_X)
            
            p_handle = patch(patch_X{p}, patch_Y{p}, [.8 .8 .8], 'FaceColor', [.8 .8 .8], 'EdgeColor', 'none');
            
            alpha(p_handle, .625)
            
        end
        
    end
    
    plot([t_interval(1) t_interval(end)], [30 30], ':w', 'LineWidth', 2)
    
    plot([t_interval(1) t_interval(end)], [8 8], ':w', 'LineWidth', 2)
    
    for c = 1:3
    
        plot(repmat(trial_times{c}, 2, 1), repmat([0; 80], 1, size(trial_times{c}, 2)), time_colors{c}, 'LineWidth', 1.5)
        
    end
    
    subplot(6, 2, 8 + ch) % subplot(4, 1, 3) % subplot(4, 2, 6 + pd)
    
    plot(t_interval, PD_dec(t_index, chan_index), 'k')
    
    hold on
    
    spikes_w_nans = Spike_indicator(t_index, chan_index);
    
    spikes_w_nans(spikes_w_nans == 0) = nan;
    
    plot(t_interval, PD_dec(t_index, chan_index).*spikes_w_nans, 'vk')
    
    set(gca, 'FontSize', 16)
    
    box off
    
    axis tight
    
    % trial_five_lims(ch, :) = ylim;
    
    for c = 1:3
    
        plot(repmat(trial_times{c}, 2, 1), repmat(ylim', 1, size(trial_times{c}, 2)), time_LFP_colors{c}, 'LineWidth', 1.5)
        
    end
    
    ylabel('LFP (mV)', 'FontSize', 16)
    
    %% Plotting LFP for 2 seconds.
    
    if isempty(zoom_loc), zoom_loc = mean(trial_borders); end
    
    sub_start_time = zoom_loc - 1/60; sub_end_time = zoom_loc + 1/60;
    
    t_sub_index = t >= sub_start_time & t <= sub_end_time;
    
    t_sub_interval = t(t_sub_index)*60;
    
    for c = 1:3
    
        subtimes{c} = times{c}(times{c} >= sub_start_time & times{c} <= sub_end_time)*60;
        
    end
    
    subplot(6, 2, 10 + ch)
    
    plot(t_sub_interval, PD_dec(t_sub_index, chan_index), 'k')
    
    set(gca, 'FontSize', 16)
    
    hold on
    
    sub_spikes_w_nans = Spike_indicator(t_sub_index, chan_index);
    
    sub_spikes_w_nans(sub_spikes_w_nans == 0) = nan;
    
    plot(t_sub_interval, PD_dec(t_sub_index, chan_index).*sub_spikes_w_nans, 'vk')
    
    set(gca, 'FontSize', 16)
    
    box off
    
    axis tight
    
    % zoom_lims(ch, :) = ylim;
    
    for c = 1:3
    
        plot(repmat(subtimes{c}, 2, 1), repmat(ylim', 1, size(subtimes{c}, 2)), time_LFP_colors{c}, 'LineWidth', 1.5)
    
    end
    
    xlabel('Time (sec.)', 'FontSize', 16)
    
    ylabel('LFP (mV)', 'FontSize', 16)
    
end

try
    
    save_as_pdf(gcf, [folder, '/', folder, '_spec_for_paper_altogether'])
    
    save_as_eps(gcf, [folder, '/', folder, '_spec_for_paper_altogether'])
    
catch
    
    save(gcf, [folder, '/', folder, '_spec_for_paper_altogether.fig'])
    
end

end
