function PD_beta_epochs_rel_infusion_roseplot_by_datapoint_individual(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

load(subject_mat)

no_subs = length(folders);

f_bins = 9:2:31; no_f_bins = length(f_bins) - 1;

f_centers = (f_bins(1:(end-1)) + f_bins(2:end))/2;

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

All_MR_mat = nan(no_f_bins, 2, no_subs, no_channels, 2);

All_conf_mat = nan(no_f_bins, 2, no_subs, no_channels, 2);

h_bins = {linspace(-pi, pi, 50), linspace(8, 32, 50)};

for i = 1:2

    h_centers{1} = (h_bins{i}(2:end) + h_bins{i}(1:(end - 1)))/2;

end
    
All_histograms = nan(50, 50, no_subs, 2, no_channels, 2);

All_slopes = nan(no_subs, 2, no_channels, 2);

All_intercepts = nan(no_subs, 2, no_channels, 2);

for ch = 1:no_channels
            
    all_beta_data = load([all_beta_name{ch},'_pbf_dp.txt']);
        
    if isempty(all_beta_data)
        
        all_beta_data = nan(1, 6);
        
    end
    
    all_subjects = all_beta_data(:, end);
        
    for s = 1:no_subs
        
        subject = str2double(folders{s});
        
        sub_indicator = all_subjects == subject;
        
        if ~isempty(sub_indicator)
            
            all_pd_index = all_beta_data(sub_indicator, 1);
            
            all_Fs = all_beta_data(sub_indicator, 3:4);
            
            all_Fc = categorize_freq(all_Fs, f_bins);
            
            all_Pds = all_beta_data(sub_indicator, 5);
            
        else
            
            all_pd_index = nan; all_Fs = nan(1, 2); all_Fc = nan(1, 2); all_Pds = nan; all_subjects = nan;
            
        end
        
        for ch1 = 1:2
            
            %% Computing Slope of Dphi/F.
            
            for pd = 1:length(pd_label)
                
                if ~isempty(all_Pds(all_pd_index == pd))
               
                    % [line_params, ~] = polyfit(all_Fs(all_pd_index == pd, ch1), all_Pds(all_pd_index == pd), 1);
               
                    [line_params, ~] = circ_on_linear(all_Fs(all_pd_index == pd, ch1), all_Pds(all_pd_index == pd), 100, 1);
                
                else
                    
                    line_params = [nan nan];
                    
                end
                
                All_slopes(s, pd, ch, ch1) = line_params(1);
                
                All_intercepts(s, pd, ch, ch1) = line_params(2);
                
            end
            
            % [pval, stats] = circ_hktest(all_Pds, all_pd_index, all_Fc, 1, {period_label{:}, f_labels{:}});
            
            figure(index)
            
            %% Plotting 2d histogram by period (pre- vs. post-infusion).
            
            for pd = 1:length(pd_label)
                
                figure(index)
                
                subplot(4, 2, pd)
                
                if ~isempty(all_Pds(all_pd_index == pd))
                    
                    [histogram, bins] = hist3([all_Pds(all_pd_index == pd) all_Fs(all_pd_index == pd, ch1)], 'Edges', h_bins);
                    
                else
                    
                    histogram = nan(50, 50); bins{1} = nan(1, 50); bins{2} = nan(1, 50);
                    
                end
                
                All_histograms(:, :, s, pd, ch, ch1) = histogram;
                
                imagesc(bins{2}, [bins{1} (bins{1} + 2*pi)], repmat(histogram, 2, 1))
                
                axis xy
                
                xlim([10 30])
                
                xlabel('Frequency (Hz)')
                
                ylabel('Phase Lag (rad)')
                
                title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];[num2str(subject), ' Phase Lag by ', chan_labels{ch1}, ' Freq.']})
                
                freezeColors
                
            end
            
            %% Plotting rose plots by period (pre- vs. post-infusion).
            
            MR_mat = nan(no_f_bins, 2); conf_mat = nan(no_f_bins, 2); no_dps = nan(no_f_bins, 2);
            
            for pd = 1:length(pd_label)
                
                figure(index)
                
                subplot(4, 2, 2 + pd)
                
                [MR_mat(:, pd), ~, ~, conf_mat(:, pd)] = rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
                
                title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];[num2str(subject), ' Phase Lag by ', chan_labels{ch1}, ' Freq.']})
                
                % figure(1)
                %
                % subplot(4, 4, (ch-1)*(2 + length(pd_label)) + (pd-1)*2 + ch1)
                %
                % rose_plot(all_Pds(all_pd_index == pd), all_Fs(all_pd_index == pd, ch1), 20, f_bins);
                %
                % title({[chan_labels{ch}, ' High Beta Blocks, ', period_label{pd}];['Phase Lag by ', chan_labels{ch1}, ' Freq.']})
                
            end
            
            All_MR_mat(:, :, s, ch, ch1) = MR_mat; All_conf_mat(:, :, s, ch, ch1) = conf_mat;
            
            %% Testing phases pre- vs. post-infusion.
            
            conc_pval = nan(no_f_bins, 1); angle_pval = nan(no_f_bins, 1);
            
            for f = 1:no_f_bins
                
                if ~isnan(any(MR_mat(f, :)))
                    
                    phi_pre = all_Pds(all_pd_index == 1 & all_Fc(:, ch1) == f);
                    
                    phi_post = all_Pds(all_pd_index == 2 & all_Fc(:, ch1) == f);
                    
                    no_dps(f, 1) = size(phi_pre, 1); no_dps(f, 2) = size(phi_post, 1);
                    
                    conc_pval(f) = circ_ktest(phi_pre, phi_post);
                    
                    %             if conc_pval(f) > 0.01/(length(f_bins) - 1)
                    
                    if ~isempty(phi_pre) && ~isempty(phi_post)
                        
                        angle_pval(f) = circ_wwtest(phi_pre, phi_post);
                        
                    else
                        
                        angle_pval(f) = 1;
                        
                    end
                    
                    %             else
                    %
                    %                 angle_pval(f) = circ_cmtest(phi_pre, phi_post);
                    %
                    %             end
                    
                else
                    
                    conc_pval(f) = 1;
                    
                    angle_pval(f) = 1;
                    
                end
                
            end
            
            % Bonferroni correcting p-values.
            conc_pval = min(conc_pval*no_f_bins, 1); angle_pval = min(angle_pval*no_f_bins, 1);
            
            %% Testing phases of frequency pairs.
            
            f_pairs = nchoosek(1:no_f_bins, 2);
            
            no_f_pairs = size(f_pairs, 1);
            
            f_conc_pval = nan(no_f_pairs, 2); f_angle_pval = nan(no_f_pairs, 2);
            
            for fp = 1:no_f_pairs
                
                for pd = 1:2
                    
                    if ~isnan(MR_mat(f_pairs(fp, 1), pd)) && ~isnan(MR_mat(f_pairs(fp, 2), pd))
                        
                        phi1 = all_Pds(all_pd_index == pd & all_Fc(:, ch1) == f_pairs(fp, 1));
                        
                        phi2 = all_Pds(all_pd_index == pd & all_Fc(:, ch1) == f_pairs(fp, 2));
                        
                        f_conc_pval(fp, pd) = circ_ktest(phi1, phi2);
                        
                        if ~isempty(phi1) && ~isempty(phi2)
                            
                            f_angle_pval(fp, pd) = circ_wwtest(phi1, phi2);
                            
                        else
                            
                            f_angle_pval(fp, pd) = 1;
                            
                        end
                        
                    else
                        
                        f_conc_pval(fp, pd) = 1;
                        
                        f_angle_pval(fp, pd) = 1;
                        
                    end
                    
                end
                
            end
            
            % Bonferroni correcting tests.
            f_conc_pval = min(f_conc_pval*no_f_pairs, 1); f_angle_pval = min(f_angle_pval*no_f_pairs, 1);
            
            figure(index)
            
            % %% Testing phases against zero.
            % 
            % zero_test = nan(no_f_bins, 2);
            % 
            % for f = 1:no_f_bins
            % 
            %     for pd = 1:2
            % 
            %         phi = all_Pds(all_pd_index == pd & all_Fc(:, ch1) == f);
            % 
            %         if ~isempty(phi)
            % 
            %             zero_test(f, pd) = circ_mtest(phi, 0);
            % 
            %         end
            % 
            %     end
            % 
            % end
            % 
            % % Bonferroni correcting p-values.
            % zero_test = min(zero_test*2*no_f_bins, 1);
            
            %% Plotting number of datapoints pre vs. post by freq.
            
            subplot(4, 3, 2*3 + 1)
            
            colormap('summer')
            
            bar(no_dps)
            
            title('Datapoints by Freq.')
            
            % legend(period_label, 'Location', 'NorthEast')
            
            set(gca, 'XTickLabel', f_centers)
            
            freezeColors
            
            %% Plotting concentration pre vs. post by freq.
            
            subplot(4, 3, 2*3 + 2)
            
            colormap('summer')
            
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
            
            freezeColors
            
            %% Plotting phase angle pre vs. post by freq.
            
            subplot(4, 3, 2*3 + 3)
            
            colormap('summer')
            
            h = barwitherr(conf_mat, angle(MR_mat));
            
            bar_pos = get_bar_pos(h);
            
            bar_pairs = {};
            
            for f = 1:no_f_bins
                
                if angle_pval(f) < 0.05
                    
                    bar_pairs = {bar_pairs{:}, [bar_pos(2*f - 1), bar_pos(2*f)]};
                    
                end
                
            end
            
            sigstar(bar_pairs, angle_pval(angle_pval < 0.05))
            
            title([num2str(subject),' Mean Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'])
            
            % legend(period_label, 'Location', 'SouthEast')
            
            set(gca, 'XTickLabel', f_centers)
            
            freezeColors
            
            %% Plotting concentration by frequency, pre and post.
            
            c_order = [linspace(1,0,no_f_bins); abs(linspace(1,0,no_f_bins)-.5); linspace(0,1,no_f_bins)]';
            
            colormap(c_order)
            
            subplot(4, 2, 3*2 + 1)%3, 2*3 + 1)
            
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
            
            % cbfreeze(h)
            
            set(gca, 'XTickLabel', period_label)
            
            %% Plotting phase angle by frequency, pre and post.
            
            subplot(4, 2, 3*2 + 2)%3, 2*3 + 2)
            
            h = barwitherr(conf_mat', angle(MR_mat)');
            
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
            
            title([num2str(subject),' Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'])
            
            colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
            
            set(gca, 'XTickLabel', period_label)
            
            % %% Plotting phase angle by frequency, pre and post, only if significantly different from zero.
            % 
            % subplot(4, 3, 2*3 + 3)
            % 
            % real_angles = angle(MR_mat);
            % 
            % real_angles(zero_test == 0) = nan;
            % 
            % h = bar(real_angles');
            % 
            % title({[num2str(subject), ' Phase Angle (', chan_labels{1}, ' - ', chan_labels{2}, ') by Freq.'];'Significantly Different from Zero'})
            % 
            % h = colorbar('YTick',1:no_f_bins,'YTickLabel',f_labels);
            % 
            % % cbfreeze(h)
            % 
            % set(gca, 'XTickLabel', period_label)
            % 
            %%
            
            save_as_pdf(index, [folders{s},'/',prefixes{s},'_',par_name,'_',ch_label{ch},'_by_ch',num2str(ch1),'_beta_ri_rose_dp'])
            
            index = index + 1;
            
        end
        
    end
    
end

save([subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_','_beta_ri_rose_dp.mat'],'All_slopes','All_intercepts','All_MR_mat','All_conf_mat','All_histograms','bins')

close('all')

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