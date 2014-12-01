function PD_outliers_pmtm(subjects_mat, epoch_secs)

load(subjects_mat)

no_figures = zeros(length(folders), 1);

for fo = 1:length(folders)
    
    clear BP_art_t
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    basetime = basetimes(fo);
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name, '_all_channel_data_dec.mat'])

    epoch_length = epoch_secs*sampling_freq;
   
    All_data = load([subj_name, '_', num2str(epoch_secs), 's_epoch_pmtm.mat']);
    
    BP_t = All_data.t;
    
    t = (1:length(PD_dec))/sampling_freq - basetime;
    
    BP_data = getfield(All_data, 'BP');
    
    bands = getfield(All_data, 'bands');
    
    for b = 1:length(bands), band_labels{b} = sprintf('%.0f - %.0f Hz', bands(b, :)); end
    
    %% Marking segments with high total power.
    
    artifact_plot = nan(size(PD_dec));
        
    artifact_indicator = zeros(length(BP_t), 2);
    
    BP_out = nan(size(BP_data));
    
    for ch = 1:2
        
        BP_zs = zscore(BP_data(:, :, ch));
        
        BP_out(:, ch) = any(BP_zs > 5, 2);
        
    end
    
    BP_out = BP_out(:, 1) | BP_out(:, 2);
    
    BP_out_index = find(BP_out);
    
    BP_out_t = BP_t(BP_out);
    
    no_outliers = length(BP_out_t);
    
    for o = 1:no_outliers
        
        segment_start = round((BP_out_t(o) + basetime)*sampling_freq - ceil(epoch_length/2));
        
        segment_start = max(segment_start, 1);
        
        segment_end = round((BP_out_t(o) + basetime)*sampling_freq + ceil(epoch_length/2));
        
        segment_end = min(segment_end, length(PD_dec));
        
        figure(1), clf
        
        subplot(1,2,1)
        
        plot(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :))
        
        axis tight
        
        title('LFP')
        
        subplot(2,2,2)
        
        bar(BP_zs(BP_out_index(o), :))
        
        title('z-Scored Band Power')
        
        set(gca, 'XTick', 1:length(bands), 'XTickLabel', band_labels)
        
        subplot(2, 2, 4)
        
        plot(t, PD_dec)
        
        hold on
        
        plot(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :), 'r')
        
        axis tight
        
        title('LFP')
        
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
            
            artifact_indicator(BP_out_index(o), ch) = 1;
            
            artifact_plot(segment_start:segment_end, ch) = PD_dec(segment_start:segment_end, ch);
            
        elseif strcmp(button, 'Cancel')
            
            return
            
        end
        
    end
    
    artifact_indicator = sum(artifact_indicator, 2) > 0;
    
    BP_art_t = BP_t(artifact_indicator);
    
    no_arts = sum(artifact_indicator);
        
    save([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.mat'], 'BP_art_t', 'artifact_indicator', 'artifact_plot')
    
    figure
    
    plot(t, [PD_dec artifact_plot])
    
    axis tight
    
    title([folder, ', Artifacts'])
    
    save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts'])
    
    rows = 3;
    
    f_no = 1;
    
    for a = 1:no_arts
        
        segment_start = round((BP_art_t(a) + basetime)*sampling_freq - ceil(epoch_length/2));
        
        segment_start = max(segment_start, 1);
        
        segment_end = round((BP_art_t(a) + basetime)*sampling_freq + ceil(epoch_length/2));
        
        segment_end = min(segment_end, length(PD_dec));
    
        if mod(a, rows^2) == 1
    
            figure
            
            f_no = f_no + 1;
    
        end
    
        subplot(rows, rows, mod(a, rows^2) + (mod(a, rows^2) == 0)*(rows^2))
    
        plot(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :))
        
        axis tight
        
        title([folder, ', Artifact ', num2str(a)])
        
        if mod(a, rows^2) == 0
            
            save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts_', num2str(f_no)])
    
        end
        
    end
    
    no_figures(fo) = f_no;
    
end

close('all')

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
   
    open([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts.fig'])
    
    for f_no = 2:no_figures(fo)
        
        open([subj_name, '_', num2str(epoch_secs), 's_pmtm_artifacts_', num2str(f_no), '.fig'])
        
    end
    
end