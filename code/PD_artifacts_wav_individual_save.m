function PD_artifacts_wav_individual_save(folder, prefix, basetime, sd_lim)

clear BP_art_t

subj_name = [folder,'/',prefix];

outliers_name = [subj_name, '_wav_BP_', num2str(sd_lim), 'sd_outliers'];

load([subj_name, '_all_channel_data_dec.mat'])

All_data = load([subj_name, '_wt_BP.mat']);

freqs = All_data.freqs;

bands = All_data.bands;

BP_t = All_data.t;

t = (1:length(PD_dec))/sampling_freq - basetime;

BP_data = All_data.BP;

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'])
    
    [~, laser_nans] = indicator_to_nans(double(laser_transitions), sampling_freq, freqs, linspace(3, 21, 200), bands);
    
    BP_data(logical(laser_nans)) = nan;
    
end

for b = 1:length(bands), band_labels{b} = sprintf('%.0f - %.0f Hz', bands(b, :)); end

%% Marking segments with high total power.

artifact_plot = nan(size(PD_dec));

artifact_indicator = zeros(length(BP_t), 2);

BP_out = nan(size(BP_data));

for ch = 1:2
    
    BP_zs = nanzscore(BP_data(:, :, ch));
    
    BP_out(:, ch) = any(BP_zs > sd_lim, 2);
    
end

BP_out = BP_out(:, 1) | BP_out(:, 2);

BP_out_blocks = index_to_blocks(BP_out);

no_outliers = size(BP_out_blocks, 1);

for o = 1:no_outliers
    
    outlier_start = BP_out_blocks(o, 1); outlier_end = BP_out_blocks(o, 2);
    
    artifact_indicator(outlier_start:outlier_end, ch) = 1;
    
    artifact_plot(outlier_start:outlier_end, ch) = PD_dec(outlier_start:outlier_end, ch);
    
end

artifact_indicator = sum(artifact_indicator, 2) > 0;

artifact_blocks = index_to_blocks(artifact_indicator);

BP_art_t = BP_t(artifact_indicator);

no_arts = size(artifact_blocks, 1);

save([outliers_name, '.mat'], 'BP_art_t', 'artifact_indicator', 'artifact_plot')

figure

plot(t, PD_dec)

hold on
plot(t, artifact_plot, 'r')

axis tight

title([folder, ', Artifacts'])

save_as_pdf(gcf, outliers_name)

rows = 3;

f_no = 0;

for a = 1:no_arts
    
    artifact_start = artifact_blocks(a, 1); artifact_end = artifact_blocks(a, 2);
    
    segment_size = diff(artifact_blocks(a, :), [], 2) + 1;
    
    segment_middle = mean(artifact_blocks(a, :), 2);
    
    segment_halfwidth = max(segment_size, 5*sampling_freq);
    
    segment_start = max(round(segment_middle - segment_halfwidth + 1), 1);
    
    segment_end = min(round(segment_middle + segment_halfwidth), length(PD_dec));
    
    if mod(a, rows^2) == 1
        
        figure
        
        f_no = f_no + 1;
        
    end
    
    subplot(rows, rows, mod(a, rows^2) + (mod(a, rows^2) == 0)*(rows^2))
    
    plot(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :))
    
    axis tight
    
    hold on
    
    plot(t(artifact_start:artifact_end), PD_dec(artifact_start:artifact_end, :), 'r')
    
    title([folder, ', Artifact ', num2str(a)])
    
    if mod(a, rows^2) == 0 || a == no_arts
        
        save_as_pdf(gcf, [outliers_name, '_', num2str(f_no)])
        
    end
    
end

end