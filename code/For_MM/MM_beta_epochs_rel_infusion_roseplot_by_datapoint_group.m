function MM_beta_epochs_rel_infusion_roseplot_by_datapoint_group(filenames, chan_labels, outlier_lim, sd_lim, win_size, smooth_size)

% SAMPLE CALL:
% MM_beta_epochs_rel_infusion_roseplot_by_datapoint_group({'file1.txt','file2.txt'},{'Striatum','Motor
% Ctx.'},7,2,333,20000)
% 
% 'filenames' is a cell of strings, which are the filenames of files 
% containing data for picking beta segments. The data should contain two
% channels, as columns.
% 'chan_labels' is a cell of strings, which are the labels of channels
% inside each file containing data.
% 'infusion_times' is a vector, with length the same as 'filenames', 
% containing the times (in datapoints) at which each data file switches
% from control to Parkinsonian behavior (or the time of carbachol
% infusion).
% 'outlier_lim' is the number of standard deviations beyond which a spike
% in the LFP is considered an outlier.
% 'sd_lim' is the number of standard deviations defining the cutoff of high
% beta power.
% 'win_size' is the minimum duration (in datapoints, so s*sampling_freq) 
% for which beta must be elevated above the cutoff, to be considered a 
% high beta segment.
% 'smooth_size' is the length of time (in datapoints, so s*sampling_freq)
% over which beta power is smoothed before applying the cutoff.

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

f_bins = 8:4:32; no_f_bins = length(f_bins) - 1;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

% f_labels = textscan(num2str(f_centers), '%s', 'delimiter', ' ');
% f_labels = cellstr(f_labels{1});
% f_labels = f_labels(1:2:end);

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

% ch_index = {1, 2, 1:2, 1, 2, 1, 2};

ch_label = {'ch1', 'ch2', 'ch1_ch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

% Creating name to call the analysis by - either name of single file in
% 'filenames', or 'All'. NB: will be over-written each time you run a
% multi-file analysis, so should be renamed after running.

if length(filenames) == 1
    
    file_label = filenames{1};
    
else

    file_label = 'All';
    
end

% Name for each channel.

all_beta_name = cell(no_channels, 1);

for ch = 1:no_channels
    
    all_beta_name{ch} = [file_label,'_',par_name,'_beta_',ch_label{ch}];
    
end

% figure(1)

index = 1;%2;

for ch = 1:no_channels
            
    all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
    
    all_pd_index = all_beta_data(:,1);
    
    all_Fs = all_beta_data(:,3:4);
    
    all_Fc = categorize_freq(all_Fs, f_bins);
    
    all_Pds = all_beta_data(:,5);
           
    for ch1 = 1:2
        
        % [pval, stats] = circ_hktest(all_Pds, all_pd_index, all_Fc, 1, {period_label{:}, f_labels{:}});
        
        figure(index)
        
        MR_mat = nan(no_f_bins, 2); no_dps = nan(no_f_bins, 2);
        
        for pd = 1:length(pd_label)
            
            figure(index)
            
            subplot(3, 2, pd)
            
            MR_mat(:, pd) = rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
            
            title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Phase Lag by ', chan_labels{ch1}, ' Freq.']})
            
            % figure(1)
            % 
            % subplot(3, 4, (ch-1)*(2 + length(pd_label)) + (pd-1)*2 + ch1)
            % 
            % rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
            % 
            % title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Phase Lag by ', chan_labels{ch1}, ' Freq.']})
                
        end
        
        conc_pval = nan(no_f_bins, 1); angle_pval = nan(no_f_bins, 1);
        
        for f = 1:no_f_bins
            
            phi_pre = all_Pds(all_pd_index == 1 & all_Fc(:, ch1) == f);
            
            phi_post = all_Pds(all_pd_index == 2 & all_Fc(:, ch1) == f);
            
            no_dps(f, 1) = size(phi_pre, 1); no_dps(f, 2) = size(phi_post, 1);
            
            conc_pval(f) = circ_ktest(phi_pre, phi_post);
            
%             if conc_pval(f) > 0.01/(length(f_bins) - 1)
            
            angle_pval(f) = circ_wwtest(phi_pre, phi_post);
                
%             else
%                 
%                 angle_pval(f) = circ_cmtest(phi_pre, phi_post); 
%                 
%             end
            
        end
        
        conc_pval = min(conc_pval*no_f_bins, 1); angle_pval = min(angle_pval*no_f_bins, 1);
        
        f_pairs = nchoosek(1:no_f_bins, 2);
        
        no_f_pairs = size(f_pairs, 1);
        
        f_conc_pval = nan(no_f_pairs, 2); f_angle_pval = nan(no_f_pairs, 2);
        
        for fp = 1:no_f_pairs
            
            for pd = 1:2
                
                phi1 = all_Pds(all_pd_index == pd & all_Fc(:, ch1) == f_pairs(fp, 1));
                
                phi2 = all_Pds(all_pd_index == pd & all_Fc(:, ch1) == f_pairs(fp, 2));
                
                f_conc_pval(fp, pd) = circ_ktest(phi1, phi2);
                
                f_angle_pval(fp, pd) = circ_wwtest(phi1, phi2);
                
            end
        
        end
        
        f_conc_pval = min(f_conc_pval*no_f_pairs, 1); f_angle_pval = min(f_angle_pval*no_f_pairs, 1);
            
        figure(index)
        
        %% Plotting number of datapoints pre vs. post by freq.
        
        subplot(3, 3, 3 + 1)
        
        bar(no_dps)
        
        title('Datapoints by Freq.')
        
        legend(period_label, 'Location', 'NorthEast')
        
        set(gca, 'XTickLabel', f_centers)
        
        %% Plotting concentration pre vs. post by freq.
        
        subplot(3, 3, 3 + 2)
        
        h = bar(abs(MR_mat));
        
        bar_pos = get_bar_pos(h);
        
        bar_pairs = {};
        
        for f = 1:no_f_bins
           
            if conc_pval(f) < 0.05
               
                bar_pairs = {bar_pairs{:}, [bar_pos(2*f - 1), bar_pos(2*f)]};
                
            end
            
        end
           
        sigstar(bar_pairs, conc_pval(conc_pval < 0.05))
        
        title('Phase Concentration by Freq.')
        
        % legend(period_label, 'Location', 'NorthWest')
        
        set(gca, 'XTickLabel', f_centers)
        
        %% Plotting phas angle pre vs. post by freq.
        
        subplot(3, 3, 3 + 3)
        
        h = bar(angle(MR_mat));
        
        bar_pos = get_bar_pos(h);
        
        bar_pairs = {};
        
        for f = 1:no_f_bins
           
            if angle_pval(f) < 0.05
               
                bar_pairs = {bar_pairs{:}, [bar_pos(2*f - 1), bar_pos(2*f)]};
                
            end
            
        end
        
        sigstar(bar_pairs, angle_pval(angle_pval < 0.05))
        
        title('Mean Phase Angle (Striatum - Motor) by Freq.')
        
        % legend(period_label, 'Location', 'SouthEast')
        
        set(gca, 'XTickLabel', f_centers)
        
        %% Plotting concentration by frequency, pre and post.
        
        subplot(3, 2, 4 + 1)
        
        h = bar(abs(MR_mat)');
        
        bar_pos = get_bar_pos(h);
        
        f_bar_pairs = {};
        
        f_conc_indicator = nan(size(f_conc_pval));
        
        for pd = 1:2
            
            % Choose whichever is smaller - significant pairs or
            % insignificant pairs.
            if sum(f_conc_pval(:, pd) < 0.05) <= no_f_pairs/2
                
                f_conc_indicator(:, pd) = f_conc_pval(:, pd) < 0.05;
                
            else
                
                f_conc_indicator(:, pd) = f_conc_pval(:, pd) >= 0.05;
                
            end
            
            for fp = 1:no_f_pairs
               
                if f_conc_indicator(fp, pd)
                    
                    f_bar_pairs = {f_bar_pairs{:}, [bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 1)), bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 2))]};
                    
                end
                
            end
            
        end
        
        f_conc_pval = reshape(f_conc_pval, 2*no_f_pairs, 1);
        
        f_conc_indicator = reshape(f_conc_indicator, 2*no_f_pairs, 1);
        
        sigstar(f_bar_pairs, f_conc_pval(f_conc_indicator == 1)')
        
        title('Phase Concentration by Freq.')
        
        h = colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);

        cbfreeze(h)
        
        set(gca, 'XTickLabel', period_label)
        
        %% Plotting phase angle by frequency, pre and post.
        
        subplot(3, 2, 4 + 2)
        
        h = bar(angle(MR_mat)');
        
        bar_pos = get_bar_pos(h);
        
        f_bar_pairs = {};
        
        f_angle_indicator = nan(size(f_angle_pval));
        
        for pd = 1:2
            
            % Choose whichever is smaller - significant pairs or
            % insignificant pairs.
            if sum(f_angle_pval(:, pd) < 0.05) <= no_f_pairs/2
                
                f_angle_indicator(:, pd) = f_angle_pval(:, pd) < 0.05;
                
            else
                
                f_angle_indicator(:, pd) = f_angle_pval(:, pd) >= 0.05;
                
            end
            
            for fp = 1:no_f_pairs
                
                if f_angle_indicator(fp, pd)
                    
                    f_bar_pairs = {f_bar_pairs{:}, [bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 1)), bar_pos((pd - 1)*no_f_bins + f_pairs(fp, 2))]};
                    
                end
                
            end
            
        end
        
        f_angle_pval = reshape(f_angle_pval, 2*no_f_pairs, 1);
        
        f_angle_indicator = reshape(f_angle_indicator, 2*no_f_pairs, 1);
        
        sigstar(f_bar_pairs, f_angle_pval(f_angle_indicator == 1)')
        
        title('Phase Angle (Striatum - Motor) by Freq.')
        
        h = colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);

        cbfreeze(h)
        
        set(gca, 'XTickLabel', period_label)
        
        %%
        
        save_as_pdf(index, [filenames(1:(end-length('_subjects.mat'))),'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_ri_rose_dp'])
        
        index = index + 1;
        
    end
    
end

save_as_pdf(gcf,[filenames(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_ri_rose_dp'])

end

function F_c = categorize_freq(F, f_bins)
    
    % Replaces continuous instantaneous frequency data with categorical frequency data.

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

    % Getting x-axis position of bars from barplot.

    for i = 1:length(handle)

        x = get(get(handle(i), 'children'), 'xdata');

        x = mean(x([1 3],:));

        pos_bars(i,:) = x;

    end

end