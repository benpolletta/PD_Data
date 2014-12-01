function PD_artifacts_wav(subjects_mat, sd_lim)

load(subjects_mat)

no_figures = zeros(length(folders), 1);

for fo = 1:length(folders)
    
    clear BP_art_t
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'])
   
    All_data = load([subj_name, '_wt_BP.mat']);
    
    BP_t = All_data.t;
    
    t = (1:length(PD_dec))/sampling_freq - basetime;
    
    BP_data = All_data.BP;
    
    bands = All_data.bands;
    
    for b = 1:length(bands), band_labels{b} = sprintf('%.0f - %.0f Hz', bands(b, :)); end
    
    %% Marking segments with high total power.
    
    artifact_plot = nan(size(PD_dec));
        
    artifact_indicator = zeros(length(BP_t), 2);
    
    BP_out = nan(size(BP_data));
    
    for ch = 1:2
        
        BP_zs = zscore(BP_data(:, :, ch));
        
        BP_out(:, ch) = any(BP_zs > sd_lim, 2);
        
    end
    
    BP_out = BP_out(:, 1) | BP_out(:, 2);
    
    BP_out_blocks = index_to_blocks(BP_out);
    
    no_outliers = size(BP_out_blocks, 1);
    
    for o = 1:no_outliers
        
        outlier_start = BP_out_blocks(o, 1); outlier_end = BP_out_blocks(o, 2);
        
        segment_size = diff(BP_out_blocks(o, :), [], 2) + 1;
        
        segment_middle = mean(BP_out_blocks(o, :), 2);
        
        segment_halfwidth = max(segment_size, 1*sampling_freq);
        
        segment_start = max(round(segment_middle - segment_halfwidth + 1), 1);
        
        segment_end = min(round(segment_middle + segment_halfwidth), length(PD_dec));
        
        figure(1), clf
        
        subplot(2, 1, 1)
        
        [ax, ~, h2] = plotyy(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :), t(segment_start:segment_end), BP_data(segment_start:segment_end, :, ch));
        
        box off
        
        hold(ax(1), 'on'), hold(ax(2), 'on')
        
        plot(ax(1), t(outlier_start:outlier_end), PD_dec(outlier_start:outlier_end, :), 'r', 'LineWidth', 2)
        
        axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
        
        % set(ax,{'ycolor'},{'black';'black'})
        
        title([folder, ', Outlier ', num2str(o)])
        
        set(get(ax(1),'YLabel'), 'String', 'LFP')
        
        set(get(ax(2),'YLabel'), 'String', 'Band Power')
        
        legend(h2, band_labels)
        
        subplot(2, 1, 2)
        
        plot(t, PD_dec)
        
        hold on
        
        plot(t(outlier_start:outlier_end), PD_dec(outlier_start:outlier_end, :), 'r')
        
        axis tight
        
        title('LFP (Entire Recording)')
        
        % [ax, ~, h2] = plotyy(t(segment_start:segment_end), PD_dec(segment_start:segment_end, ch), t(segment_start:segment_end), BP_zs(segment_start:segment_end, :));
        %
        % box off
        %
        % hold(ax(1), 'on'), hold(ax(2), 'on')
        %
        % BP_multiplier = BP_out(segment_start:segment_end, :);
        %
        % BP_multiplier(BP_multiplier == 0) = nan;
        %
        % plot(ax(2), t(segment_start:segment_end), BP_multiplier.*BP_zs(segment_start:segment_end, :), 'LineWidth', 2)
        %
        % axis(ax(1), 'tight'), axis(ax(2), 'tight'), xlim(ax(2),xlim(ax(1)))
        %
        % set(ax,{'ycolor'},{'black';'black'})
        %
        % set(get(ax(1),'YLabel'), 'String', 'LFP')
        %
        % set(get(ax(2),'YLabel'), 'String', 'Band Power')
        %
        % legend(h2, band_labels, 'Location', 'SouthWest')
        
        button = questdlg('Is this an artifact?');
        
        if strcmp(button, 'Yes')
            
            artifact_indicator(outlier_start:outlier_end, ch) = 1;
            
            artifact_plot(outlier_start:outlier_end, ch) = PD_dec(outlier_start:outlier_end, ch);
            
        elseif strcmp(button, 'Cancel')
            
            return
            
        end
        
    end
    
    artifact_indicator = sum(artifact_indicator, 2) > 0;
    
    artifact_blocks = index_to_blocks(artifact_indicator);
    
    BP_art_t = BP_t(artifact_indicator);
    
    no_arts = size(artifact_blocks, 1);
        
    save([subj_name, '_wav_BP_', num2str(sd_lim), 'sd_artifacts.mat'], 'BP_art_t', 'artifact_indicator', 'artifact_plot')
    
    figure
    
    plot(t, [PD_dec artifact_plot])
    
    axis tight
    
    title([folder, ', Artifacts'])
    
    save_as_pdf(gcf, [subj_name, '_wav_BP_', num2str(sd_lim), 'sd_artifacts'])
    
    rows = 3;
    
    f_no = 1;
    
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
        
        if mod(a, rows^2) == 0
            
            save_as_pdf(gcf, [subj_name, '_wav_BP_', num2str(sd_lim), 'sd_artifacts_', num2str(f_no)])
    
        end
        
    end
    
    no_figures(fo) = f_no;
    
end

close('all')

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    open([subj_name, '_wav_BP_', num2str(sd_lim), 'sd_artifacts.fig'])
    
    for f_no = 1:no_figures(fo)
        
        open([subj_name, '_wav_BP_', num2str(sd_lim), 'sd_artifacts_', num2str(f_no), '.fig'])
        
    end
    
end