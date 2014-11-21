function PD_outliers_pmtm(subjects_mat, epoch_secs)

load(subjects_mat)

for fo = 1:length(folders)
    
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
    
    %% Marking segments with high total power.
    
    outliers_plot = nan(size(PD_dec));
    
    for ch = 1:2
        
        BP_zs = zscore(BP_data(:, :, ch));
        
        BP_out = any(BP_zs > 5, 2);
        
        BP_out_t = BP_t(BP_out);
        
        no_outliers = length(BP_out_t);
        
        for o = 1:no_outliers
            
            segment_start = round((BP_out_t(o) + basetime)*sampling_freq - ceil(epoch_length/2));
            
            segment_start = max(segment_start, 1);
            
            segment_end = round((BP_out_t(o) + basetime)*sampling_freq + ceil(epoch_length/2));
            
            segment_end = min(segment_end, length(PD_dec));
            
            outliers_plot(segment_start:segment_end, ch) = PD_dec(segment_start:segment_end, ch);
            
        end
        
    end
        
    save([subj_name, '_', num2str(epoch_secs), 's_pmtm_outliers.mat'], 'BP_out_t')
    
    figure
    
    plot(t, [PD_dec outliers_plot])
    
    axis tight
    
    title([folder, ', Outliers'])
    
    save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_pmtm_outliers_1'])
    
    rows = 3;
    
    f_no = 1;
    
    for o = 1:no_outliers
        
        segment_start = round((BP_out_t(o) + basetime)*sampling_freq - ceil(epoch_length/2));
        
        segment_start = max(segment_start, 1);
        
        segment_end = round((BP_out_t(o) + basetime)*sampling_freq + ceil(epoch_length/2));
        
        segment_end = min(segment_end, length(PD_dec));
    
        if mod(o, rows^2) == 1
    
            figure
            
            f_no = f_no + 1;
    
        end
    
        subplot(rows, rows, mod(o, rows^2) + (mod(o, rows^2) == 0)*(rows^2))
    
        plot(t(segment_start:segment_end), PD_dec(segment_start:segment_end, :))
        
        axis tight
        
        title([folder, ', Outlier ', num2str(o)])
        
        if mod(o, rows^2) == 0
            
            save_as_pdf(gcf, [subj_name, '_', num2str(epoch_secs), 's_pmtm_outliers_', num2str(f_no)])
    
        end
        
    end
    
end