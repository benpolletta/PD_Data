function PD_laser_times_wav(subjects_mat)

load(subjects_mat)

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    data_length = length(PD_dec);
    
    laser_transitions = zeros(data_length, 1);
    
    laser_periods = zeros(data_length, 3);
    
    laser_on = laseron(fo);
    
    lasertime = lasertimes(fo);
    
    time = (1:data_length)/sampling_freq; % - laser_on;
    
    trial_length = triallength(fo);
    
    no_trials = ceil((data_length/sampling_freq - laser_on)/trial_length);
    
    for t = 1:no_trials
       
        trial_laser_on = single((laser_on + (t - 1)*trial_length)*sampling_freq);
        
        trial_laser_off = single((laser_on + (t - 1)*trial_length + lasertime)*sampling_freq);
        
        trial_off = single((laser_on + t*trial_length)*sampling_freq);
        
        laser_transitions(trial_laser_on) = 1;
        
        laser_transitions(trial_laser_off) = 1;
        
        laser_periods(trial_laser_on:trial_laser_off, 1) = 1;
        
        laser_periods(trial_laser_off + (1:(5*sampling_freq)), 2) = 1;
        
        laser_periods(trial_off - (1:(5*sampling_freq)), 3) = 1;
        
    end
    
    save([subj_name, '_wav_laser_times.mat'], 'laser_transitions', 'laser_periods')
    
    figure
    
    time_limit = inf; %120;
    
    plot_periods = laser_periods;
    
    plot_periods(plot_periods == 0) = nan;
    
    plot_periods(plot_periods == 1) = 0;
    
    plot(time(time < time_limit), [PD_dec(time < time_limit, :) plot_periods(time < time_limit, :)]) % norm_by_mat(laser_periods, PD_dec)])
    
    hold on
    
    transition_times = time(logical(laser_transitions));
    
    transition_times(transition_times > time_limit) = [];
    
    plot(repmat(transition_times, 2, 1), diag([all_dimensions(@nanmin, PD_dec) all_dimensions(@max, PD_dec)])*ones(2, length(transition_times)), 'k', 'LineWidth', 2)
    
    axis tight
    
    title([folder, ', Laser On & Off Times'])
    
    legend({chan_labels{:}, 'Laser On', 'Post-Laser', 'Pre-Laser', 'Laser Transitions'})
    
    xlabel('Time (s)')
    
    save_as_pdf(gcf, [subj_name, '_wav_laser_times'])
    
end
