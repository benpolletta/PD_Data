function PD_rel_infusion_collect_freq(subject_mat)

load(subject_mat)

sampling_freq = 1000;

smooth_winsize = 50;

period_label = {'Pre-Infusion','Post-Infusion'};

no_periods = length(period_label);

all_freq_name = [subject_mat(1:(end-length('_subject.mat'))),'freq.txt'];

if isempty(dir(all_freq_name))
    
    fid_all = fopen(all_freq_name, 'w');
    
    for fo = 1:length(folders)
        
        folder = folders{fo};
        
        prefix = prefixes{fo};
        
        base_index = basetimes(fo)*sampling_freq;
        
        subj_name = [folder,'/',prefix];
        
        load([subj_name,'_all_channel_data_dec_HAP.mat'])
        
        periods = [1 base_index; (base_index + 1) size(P,1)];
        
        no_periods = size(periods,1);
        
        freq_name = [subj_name,'_freq'];
        
        figure;
        
        h_all = nan(50, 2*no_periods, 2);
        
        b_all = nan(50, 2*no_periods, 2);
        
        for pd = 1:no_periods
            
            P_pd = unwrap(P(periods(pd, 1):periods(pd, 2), :, 3));
            
            F = diff(P_pd)/(2*pi*(1/sampling_freq));
            
            F_smooth = nan(size(P_pd,1),2);
            
            for ch1 = 1:2
                
                F_flipped = [flipud(F(1:smooth_winsize,ch1)); F(:,ch1); flipud(F((end-smooth_winsize+1):end,ch1))];
                
                % F_conv = conv(F_flipped,hann(smooth_winsize)/sum(hann(smooth_winsize)),'same');
                F_conv = conv(F_flipped, ones(smooth_winsize,1)/smooth_winsize, 'same');
                
                F_smooth(:,ch1) = F_conv(smooth_winsize + (1:size(P_pd,1)));
                
            end
            
            F_rk = tiedrank(F_smooth);
            
            fprintf(fid_all, '%f\t%f\t%f\t%f\t%f\n', [pd*ones(size(F_rk), 1) F_smooth F_rk*diag(1./max(F_rk))]');
            
            subplot(2, no_periods, pd)
            
            r = nancorr(F_smooth(:, 1), F_smooth(:, 2));
            
            [histogram, centers] = hist3(F_smooth, [50 50]);
            
            colormap('default')
            
            imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
            
            h = colorbar;
            
            title({[folder,' ',period_label{pd}];['Correlation = ',num2str(r)]})
            
            ylabel(h, 'Prop. Observed')
            
            axis xy
            
            xlabel([chan_labels{1},' Frequency'])
            
            ylabel([chan_labels{2}, ' Frequency'])
            
            for ch1 = 1:2
                
                top_q_cutoff = quantile(F_smooth(:, ch1), 0.75);
                
                top_q_indicator = F_smooth(:, ch1) >= top_q_cutoff;
                
                other_amp = F_smooth(top_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, pd, ch1) = h/sum(h);
                
                b_all(:, pd, ch1) = b;
                
                q_legend{pd} = [period_label{pd}, ', Top 25%'];
                
                bot_q_cutoff = quantile(F_smooth(:, ch1), 0.25);
                
                bot_q_indicator = F_smooth(:, ch1) <= bot_q_cutoff;
                
                other_amp = F_smooth(bot_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, 2 + pd, ch1) = h/sum(h);
                
                b_all(:, 2 + pd, ch1) = b;
                
                q_legend{2 + pd} = [period_label{pd}, ', Bottom 25%'];
                
            end
            
        end
        
        for ch1 = 1:2
            
            subplot(2, no_periods, no_periods + ch1)
            
            plot(b_all(:, :, ch1), h_all(:, :, ch1))
            
            legend(q_legend)
            
            title([chan_labels{3 - ch1}, ' Freq. Dist. for Top/Bottom Quartile of ', chan_labels{ch1}, ' Freq.'])
            
            ylabel('Proportion of Observations')
            
            xlabel([chan_labels{3 - ch1}, ' Beta Freq.'])
            
        end
        
        save_as_pdf(gcf, freq_name)
        
    end
    
    fclose(fid_all);
    
end

all_freq_data = load(all_freq_name);

all_pds = all_freq_data(:, 1);

norm_labels = {'Freq.','Freq. (Rank)'};

no_norms = length(norm_labels);

figure

for norm = 1:no_norms
    
    for pd = 1:no_periods
        
        subplot(no_norms + 1, no_periods, (norm - 1)*no_periods + pd)
        
        F_pd = all_freq_data(all_pds == pd, (norm - 1)*no_periods + (2:3));
        
        r = nancorr(F_pd(:, 1), F_pd(:, 2));
        
        [histogram, centers] = hist3(F_pd, [50 50]);
        
        colormap('default')
        
        imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
        
        h = colorbar;
        
        title({[chan_labels{1},' by ',chan_labels{2},', ',period_label{pd}];['Correlation = ',num2str(r)]})
        
        ylabel(h, 'Prop. Observed')
        
        axis xy
        
        xlabel([chan_labels{1}, ' ', norm_labels{norm}])
        
        ylabel([chan_labels{2}, ' ', norm_labels{norm}])
        
        if norm == 1
        
            for ch1 = 1:2
                
                top_q_cutoff = quantile(F_pd(:, ch1), 0.75);
                
                top_q_indicator = F_pd(:, ch1) >= top_q_cutoff;
                
                other_amp = F_pd(top_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, pd, ch1) = h/sum(h);
                
                b_all(:, pd, ch1) = b;
                
                q_legend{pd} = [period_label{pd}, ', Top 25%'];
                
                bot_q_cutoff = quantile(F_pd(:, ch1), 0.25);
                
                bot_q_indicator = F_pd(:, ch1) <= bot_q_cutoff;
                
                other_amp = F_pd(bot_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, 2 + pd, ch1) = h/sum(h);
                
                b_all(:, 2 + pd, ch1) = b;
                
                q_legend{2 + pd} = [period_label{pd}, ', Bottom 25%'];
                
            end
            
        end
        
    end
    
    if norm == 1
        
        for ch1 = 1:2
            
            subplot(no_norms + 1, no_periods, no_norms*no_periods + ch1)
            
            plot(b_all(:, :, ch1), h_all(:, :, ch1))
            
            legend(q_legend)
            
            title([chan_labels{3 - ch1}, ' Freq. Dist. for Top/Bottom Quartile of ', chan_labels{ch1}, ' Freq.'])
            
            ylabel('Proportion of Observations')
            
            xlabel([chan_labels{3 - ch1}, ' Beta Freq.'])
            
        end
        
    end
    
end

save_as_pdf(gcf, all_freq_name(1:end-4))