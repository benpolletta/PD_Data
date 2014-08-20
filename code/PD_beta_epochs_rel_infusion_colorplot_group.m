function PD_beta_epochs_rel_infusion_colorplot_group(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

f_bins = 9:2:31; no_f_bins = length(f_bins) - 1;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';

f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
f_labels = cellstr(f_labels{1});
f_labels = f_labels(1:2:end);

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

% figure(1)

index = 1;%2;

for ch = 1:no_channels
            
    all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
    
    if ~isempty(all_beta_data)
        
        all_pd_index = all_beta_data(:,1);
        
        all_Fs = all_beta_data(:,3:4);
        
        all_Fc = categorize_freq(all_Fs, f_bins);
        
        all_Pds = all_beta_data(:,5);
        
    else
       
        all_pd_index = nan; all_Fs = nan(1, 2); all_FC = nan(1, 2); all_Pds = nan;
        
    end
           
    for ch1 = 1:2
        
        % [pval, stats] = circ_hktest(all_Pds, all_pd_index, all_Fc, 1, {period_label{:}, f_labels{:}});
        
        figure(index)
        
        MR_mat = nan(no_f_bins, 2); conf_mat = nan(no_f_bins, 2); no_dps = nan(no_f_bins, 2);
            
        %% Plotting 2d histogram by period (pre- vs. post-infusion).
        
        for pd = 1:length(pd_label)
            
            figure(index)
            
            subplot(4, 2, pd)
            
            if ~isempty(all_Pds(all_pd_index == pd))
            
                [histogram, bins] = hist3([all_Pds(all_pd_index == pd) all_Fs(all_pd_index == pd, ch1)], [50 50]);
            
            else
               
                histogram = nan(50, 50); bins{1} = nan(1, 50); bins{2} = nan(1, 50);
                
            end
                
            imagesc(bins{2}, [bins{1} (bins{1} + 2*pi)], repmat(histogram, 2, 1)) %imagesc(bins{2}, bins{1}, histogram)
            
            axis xy
            
            xlim([10 30])
            
            xlabel('Frequency (Hz)')
            
            ylabel('Phase Lag (rad)')
            
            title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];[' Phase Lag by ', chan_labels{ch1}, ' Freq.']})
            
            freezeColors
            
        end
        
        save_as_pdf(index, [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_ri_rose_dp'])
        
        index = index + 1;
        
    end
    
end

save_as_pdf(gcf,[subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_ri_rose_dp'])

end

function F_c = categorize_freq(F, f_bins)
    
    [r, c] = size(F);
    
    F_c = zeros(r, c);
    
    no_f_bins = length(f_bins) - 1;
    
    for col = 1:c
        
        for f = 1:no_f_bins
            
            F_bin = F(:, col) >= f_bins(f) & F(:, col) < f_bins(f + 1);
            
            F_c(:, col) = F_c(:, col) + f*F_bin;
            
        end
        
    end
    
end

function pos_bars = get_bar_pos(handle)

for i = 1:length(handle)
    
    x = get(get(handle(i), 'children'), 'xdata');
    
    x = mean(x([1 3],:));
    
    pos_bars(i,:) = x;
    
end

end