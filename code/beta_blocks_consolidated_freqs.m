function beta_blocks_consolidated_freqs(subject_mat, peak_suffix, time_window, percent, norm, band_index, freqs, no_cycles, bands)

PD_struct = PD_initialize(subject_mat);

if isempty(freqs) && isempty(no_cycles) && isempty(bands)
    
    freqs = 1:200;
    
    no_cycles = linspace(3, 21, length(freqs));
    
    bands = [1 4; 4 8; 8 30; 30 100; 120 180; 0 200];
    
    BP_suffix = peak_suffix;
    
else
    
    BP_suffix = sprintf('_%.0f-%.0fHz_%.0f-%.0fcycles_%dbands', freqs(1), freqs(end), no_cycles(1), no_cycles(end), size(bands, 1));
    
    BP_suffix = [BP_suffix, peak_suffix];
    
end

no_bands = size(bands, 1);

[band_indices, short_band_labels, band_labels] = deal(cell(no_bands, 1));

for b = 1:no_bands
   
    band_indices{b} = freqs >= bands(b, 1) & freqs <= bands(b, 2);
    
    short_band_labels{b} = sprintf('%d-%dHz', bands(b, :));
    
    band_labels{b} = sprintf('%d - %d Hz', bands(b, :));
    
end

band_freqs = freqs(band_indices{band_index});

format = make_format(sum(band_indices{band_index}) + 2, 'f');

fid = nan(PD_struct.no_pds, PD_struct.no_chans);

for pd = 1:PD_struct.no_pds
    
    for ch = 1:PD_struct.no_chans
    
        fid(pd, ch) = fopen([PD_struct.subj_prefix, BP_suffix, '_2sd_', short_band_labels{band_index}, '_high_',...
        num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated_freqs',...
        norm, '_ch', num2str(ch), '_', PD_struct.pd_labels{pd}, '.txt'], 'w');
        
    end
    
end

for fo = 1:PD_struct.no_folders
    
    folder = PD_struct.folders{fo};
    
    prefix = PD_struct.prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    cons_name = [subj_name, BP_suffix, '_2sd_BP_high_', num2str(time_window/PD_struct.sampling_freq), 's_', num2str(percent), 'pct_consolidated'];
    
    % cons_title = [folder, ', High Power, ', num2str(percent), ' Percent of ', num2str(time_window/PD_struct.sampling_freq), ' s'];
    
    load([cons_name, '.mat'], 'BP_high_consolidated')
    
    load([subj_name, '_all_channel_data_dec.mat'])
    
    t = (1:length(PD_dec))/PD_struct.sampling_freq - PD_struct.basetimes(fo);
    
    pd_indices = get_pd_indices(t, fo, subj_name, PD_struct.subj_prefix, [], '', 3, PD_struct.no_pds);
    
    load([subj_name, '_wt.mat'], 'Spec')
    
    BP_high_cons = getfield(BP_high_consolidated, ['BP_high', PD_struct.high_type{2}]);
    
    for pd = 1:PD_struct.no_pds
        
        %% Get Spec data.
        
        interval_length = 10;
        
        BP_hc_length = sum(BP_high_cons(:, band_index, 3) & pd_indices(:, pd));
        
        if BP_hc_length > 0
            
            BP_hc_blocks = index_to_blocks(BP_high_cons(:, band_index, 3) & pd_indices(:, pd));
            
            no_blocks = size(BP_hc_blocks, 1);
            
            plot_length = BP_hc_length + (no_blocks - 2)*interval_length;
            
            Spec_hc = nan(plot_length, sum(band_indices{band_index}), 2);
            
            t_hc_plot = nan(plot_length, 1);
            
            last_index = 0;
            
            for bl = 1:no_blocks
                
                block_length = BP_hc_blocks(bl, 2) - BP_hc_blocks(bl, 1) + 1;
                
                block_indices = last_index + (bl > 1)*5 + (1:block_length);
                
                t_hc_plot(block_indices) = t(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2));
                
                Spec_hc(block_indices, :, :) = Spec(BP_hc_blocks(bl, 1):BP_hc_blocks(bl, 2), band_indices{band_index}, :);
                
                last_index = last_index + (bl > 1)*interval_length + block_length;
                
            end
            
            Freq_indices = nan(size(Spec_hc, 1), 1);
            
            [Phases_high_beta, Freqs_high_beta] = deal(nan(size(Spec_hc, 1), 2));
            
            for ch = 1:PD_struct.no_chans
                
                Spec_hc_chan = Spec_hc(:, :, ch);
                
                if strcmp(norm, '_zscore')
                    
                    Spec_zs = nanzscore(abs(Spec_hc_chan));
                    
                    [~, Freq_indices] = nanmax(Spec_zs, [], 2);
                    
                else
                    
                    [~, Freq_indices] = nanmax(abs(Spec_hc_chan), [], 2);
                    
                end
                
                Freq_indices_linear = sub2ind(size(Spec_hc_chan), (1:size(Spec_hc_chan, 1))', Freq_indices);
                
                Phases_high_beta(:, ch) = angle(Spec_hc_chan(Freq_indices_linear));
                
                Freqs_high_beta(:, ch) = band_freqs(Freq_indices);
                
                fprintf(fid(pd, ch), format, [Freqs_high_beta(:, ch) Phases_high_beta(:, ch) Spec_hc(:, :, ch)]');
                
            end
            
        else
           
            [Phases_high_beta, Freqs_high_beta] = deal(nan(0, 2));
            
            Spec_hc = nan(0, sum(band_indices{band_index}));
            
        end
        
        save([cons_name, '_', short_band_labels{band_index}, '_', PD_struct.pd_labels{pd}, '_freqs', norm, '.mat'], 'Freqs_high_beta', 'Phases_high_beta', 'Spec_hc')
        
    end
    
end

fclose('all');

end

function pd_indices = get_pd_indices(t, fo, subj_name, subj_mat_name, epoch_secs, BP_suffix, band_index, no_pds)

if ~isempty(dir([subj_name, '_wav_laser_artifacts.mat']))
    
    load([subj_name, '_wav_laser_artifacts.mat'], 'laser_periods')
    
    pd_indices = laser_periods;
    
elseif ~isempty(epoch_secs)
    
    if ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR.mat']))
        
        load([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs_by_STR.mat'])
        
        pd_indices = zeros(length(t), no_pds);
        
        for pd = 1:no_pds
            
            bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
            
            bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
            
            pd_indices(bp_max_start:bp_max_end, pd) = 1;
            
        end
        
    elseif ~isempty(dir([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat']))
        
        load([subj_mat_name, BP_suffix, '_pct_BP_high_', num2str(epoch_secs/60), '_min_secs.mat'])
        
        pd_indices = zeros(length(t), no_pds);
        
        for pd = 1:no_pds
            
            bp_max_start = All_bp_max_start(fo, 1, band_index, pd);
            
            bp_max_end = All_bp_max_end(fo, 1, band_index, pd);
            
            pd_indices(bp_max_start:bp_max_end, pd) = 1;
            
        end
        
    end
    
else
    
    if no_pds == 2
        
        pd_indices(:, 1) = t < 0;
        
        pd_indices(:, 2) = t > 0;
        
    elseif no_pds == 1
        
        pd_indices = ones(length(t), no_pds);
        
    end
    
end

end

