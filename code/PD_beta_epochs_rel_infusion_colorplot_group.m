function PD_beta_epochs_rel_infusion_colorplot_group(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

bands = [10 18; 18 22; 22 30];

no_bands = size(bands,1);

colors = [.5 .5 .5; .5 .5 .5; .5 .5 .5]; %[1 1 1; 0 0 0; .5 .5 .5; .5 .5 .5; .5 .5 .5];

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', 'Either', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

% ch_index = {1, 2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1andch2', 'ch1orch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

all_beta_name = cell(no_channels, 1);

for ch = 1:no_channels
    
    all_beta_name{ch} = [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_',ch_label{ch}];
    
end

h_bins = {linspace(-pi, pi, 50), linspace(8, 32, 50)};

for i = 1:2

    h_centers{i} = (h_bins{i}(2:end) + h_bins{i}(1:(end - 1)))/2;

end

index = 1;

for ch = 1:no_channels
            
    all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
    
    if ~isempty(all_beta_data)
        
        all_pd_index = all_beta_data(:,1);
        
        all_Fs = all_beta_data(:,3:4);
        
        all_Pds = all_beta_data(:,5);
        
    else
       
        all_pd_index = nan; all_Fs = nan(1, 2); all_Pds = nan;
        
    end
           
    for ch1 = 1:2
        
        figure(index)
        
        for pd = 1:length(pd_label)
            
            subplot(1, 2, pd)
            
            if ~isempty(all_Pds(all_pd_index == pd))
                
                [histogram, bins] = hist3([all_Pds(all_pd_index == pd) all_Fs(all_pd_index == pd, ch1)], 'Edges', h_bins);
                
            else
                
                histogram = nan(50, 50); bins{1} = nan(1, 50); bins{2} = nan(1, 50);
                
            end
            
            imagesc(bins{2}, bins{1}, histogram) %[bins{1} (bins{1} + 2*pi)], repmat(histogram, 2, 1))
            
            hold on
            
            axis xy
            
            xlim([10 30])
            
            xlabel('Frequency (Hz)')
            
            ylabel('Phase Lag (rad)')
            
            title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Phase Lag by ', chan_labels{ch1}, ' Freq.']})
            
            freezeColors
            
            %% Computing Slope of Dphi/F.
                
            line_params = nan(no_bands, 2);
            
            if ~isempty(all_Pds(all_pd_index == pd))
                
                % [line_params(1, :), ~] = polyfit(all_Fs(all_pd_index == pd, ch1), all_Pds(all_pd_index == pd), 1);
                
                % [line_params(2, :), ~] = circ_on_linear(all_Fs(all_pd_index == pd, ch1), all_Pds(all_pd_index == pd), 100, 0);
                
                for b = 1:no_bands
                    
                    freq_index = all_pd_index == pd & all_Fs(:, ch1) <= bands(b, 2) & all_Fs(:, ch1) >= bands(b, 1);
                    
                    [line_params(b, :), ~] = polyfit(all_Fs(freq_index, ch1), all_Pds(freq_index), 1); %[line_params(2 + b, :), ~]
                    
                end
                
            end
            
            % Plotting.
            
            for b = 1:no_bands %(no_bands + 2)
                
                line_x_vals = h_centers{2};
                
                % if b >= 3
                    
                line_x_vals = line_x_vals(line_x_vals <= bands(b, 2) & line_x_vals >= bands(b, 1)); %bands(b - 2, 2)
                    
                % end
                
                line_y_vals = polyval(line_params(b, :), line_x_vals);
                
                plot(line_x_vals, line_y_vals, 'Color', colors(b, :), 'LineWidth', 3)
                
            end
            
        end
        
        save_as_pdf(index, [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_cplot'])
        
        index = index + 1;
        
    end
    
end