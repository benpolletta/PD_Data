function PD_beta_blocks_rel_infusion_pre_post_PLV_plot_individual(subject_mat, peak_suffix, epoch_secs, pd_handle, band_index_for_time, band_index_for_display, freqs, no_cycles, bands, folders_left_out)

% Leave epoch_secs empty when using for optogenetics data, and enter
% '_ntrials' for the argument pd_handle.

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = ['', peak_suffix];
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end
    
% close('all')

load(subject_mat)

load('bonferroni_count')

no_folders = length(folders);

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
    
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
   
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

display_indices = band_indices{band_index_for_display};

no_chans = length(chan_labels);

no_pds = length(pd_labels);

subj_mat_name = subject_mat(1:(end - length('_subjects.mat')));

if ~isempty(folders_left_out)
    
    if ischar(folders_left_out)
    
        folder_flag = ['_no_', folders_left_out];
        
        folder_chi = ~strcmp(folders, folders_left_out);
        
    elseif iscell(folders_left_out)
        
        folder_flag = folders_left_out{1};
        
        folder_cell = folders_left_out{2}; 
        
        folder_chi = ones(1, length(folders));
        
        for fo = 1:length(folder_cell)
           
            folder_chi(strcmp(folders, folder_cell{fo})) = 0;
            
        end
        
    end
    
else
    
    folder_flag = '';
    
    folder_chi = ones(1, length(folders));
    
end

if ~isempty(epoch_secs)
    
    PLV_name = [subj_mat_name, BP_suffix, '_pct_', short_band_labels{band_index_for_time}, '_high_',...
        num2str(epoch_secs/60), '_min_secs', pd_handle, '_PLV'];
    
else
    
    PLV_name = [subj_mat_name, BP_suffix, '_pct_',...
        short_band_labels{band_index_for_time}, pd_handle, '_PLV'];
    
    epoch_secs = str2double(pd_handle(2:(end - length('trials'))));
    
end

load([PLV_name, '.mat'])

no_freqs = sum(display_indices);

array_names = {'Coh_sec', 'Coh_sec_pct', 'dP_sec'}; % , 'dP_sec_pct'};

long_array_names = {'Coherence', 'Coherence (%\Delta Baseline)', 'Phase of Coh.'}; % , 'Phase of Coh. (%\Delta Baseline)'};

no_arrays = length(array_names);

All_mean = nan(no_freqs, no_folders, no_pds, no_arrays);

%% Individual bar plots.

figure

for fo = 1:no_folders
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder, '/', prefix];
    
    for a = 1:no_arrays
        
        eval(['PLV_sec = ', array_names{a}, ';'])
        
        [PLV_mean, PLV_ci] = deal(nan(no_freqs, no_pds));
        
        for pd = 1:no_pds
            
            if strcmp(array_names{a}, 'dP_sec')
                
                PLV_mean(:, pd) = angle(nanmean(exp(sqrt(-1)*PLV_sec(display_indices, :, fo, pd)), 2));
    
                % PLV_mean(PLV_mean(:, pd) < -pi/2, pd) = 2*pi + PLV_mean(PLV_mean(:, pd) < -pi/2, pd);
                
                PLV_ci(:, pd) = circ_confmean(PLV_sec(display_indices, :, fo, pd), .05, [], [], 2)/sqrt(epoch_secs); % .05/bonferroni_count
                
            else
                
                PLV_mean(:, pd) = nanmean(PLV_sec(display_indices, :, fo, pd), 2);
                
                PLV_ci(:, pd) = norminv(1 - .05, 0, 1)*... % .05/length(freqs(display_indices)
                    nanstd(PLV_sec(display_indices, :, fo, pd), [], 2)/sqrt(epoch_secs);
                
            end
            
        end
        
        subplot(no_arrays, no_folders, (a - 1)*no_folders + fo)
        % subplot(no_folders, no_chans, (fo - 1)*no_chans + ch) % subplot(no_bands, 2, (b - 1)*2 + ch)
            
        if strcmp(array_names{a}, 'dP_sec')
            
            boundedline(freqs(display_indices), (180/pi)*PLV_mean, prep_for_boundedline((180/pi)*PLV_ci))
            
        else
            
            boundedline(freqs(display_indices), PLV_mean, prep_for_boundedline(PLV_ci))
            
        end
        
        axis tight
        
        if folder_chi(fo)
            
            All_mean(:, fo, :, a) = permute(PLV_mean, [1 3 2]);
            
            if strcmp(array_names{a}, 'dP_sec')
               
                All_mean(:, fo, 1, a + 1) = diff(PLV_mean, [], 2);
                
            end
            
        end
        
        if fo == 1
            
            title({[num2str(epoch_secs/60), ' Minutes of Densest High Power'];[folder, ', ', band_labels{b}]})
            % title({chan_labels{ch};[num2str(epoch_secs/60), ' Minutes of Densest High Power'];[folder, ', ', band_labels{b}]})
            
            ylabel(array_names{a})
            
            legend(pd_labels)
            
        else
            
            title(folder)
            % ylabel(folder)
            
        end
        
    end
    
end
    
save([PLV_name, '_data_for_plot.mat'], 'All_mean')

save_as_pdf(gcf, [PLV_name, '_', short_band_labels{band_index_for_display}, '_individual'])

%% Stats & figures treating each individual as an observation.

array_names = {'Coh_sec', 'Coh_sec_pct', 'dP_sec', 'ddP_sec'};

long_array_names = {'Coherence', 'Coherence (%\Delta Baseline)', 'Phase of Coh.', '\Delta Phase of Coh.'}; % , 'Phase of Coh. (%\Delta Baseline)'};

no_arrays = length(array_names);

figure

for a = 1:no_arrays

    [All_mean_mean, All_mean_ci] = deal(nan(length(freqs(display_indices)), no_pds));
    
    for pd = 1:no_pds
        
        if strcmp(array_names{a}((end - length('dP_sec') + 1):end), 'dP_sec')
            
            All_mean_mean(:, pd) = (180/pi)*angle(nanmean(exp(sqrt(-1)*All_mean(:, :, pd, a)), 2));
            
            All_mean_ci(:, pd) = (180/pi)*circ_confmean(All_mean(:, :, pd, a), .05, [], [], 2); % .05/length(freqs(display_indices)
            
        else
            
            All_mean_mean(:, pd) = nanmean(All_mean(:, :, pd, a), 2);
            
            All_mean_ci(:, pd) = norminv(1 - .05, 0, 1)*... % .05/length(freqs(display_indices))
                nanstd(All_mean(:, :, pd, a), [], 2)/sqrt(no_folders);
            
        end
        
    end
    
    save([PLV_name, folder_flag, '_', array_names{a}, '_data_for_plot.mat'], 'All_mean_mean', 'All_mean_ci')
    
    subplot(1, no_arrays, a)
    
    boundedline(freqs(display_indices), All_mean_mean, prep_for_boundedline(All_mean_ci))
        
    axis tight
    
    title([long_array_names{a}, ', ', num2str(epoch_secs/60), ' Minutes of Densest High Power, ', band_labels{band_index_for_time}])
    
    legend(pd_labels)
    
end

save_as_pdf(gcf, [PLV_name, folder_flag, '_', short_band_labels{band_index_for_display}, '_individual_avg'])

end
