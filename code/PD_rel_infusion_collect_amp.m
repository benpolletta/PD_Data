function PD_rel_infusion_collect_amp(subject_mat)

load(subject_mat)

sampling_freq = 1000; win_size = 333;

period_label = {'Pre-Infusion','Post-Infusion'};

no_periods = length(period_label);

all_amp_name = [subject_mat(1:(end-length('_subject.mat'))),'amp.txt'];

% if isempty(dir(all_amp_name))

fid_all = fopen(all_amp_name, 'w');

all_corr = nan(length(folders), length(-win_size:win_size), 2);

all_lags = nan(length(folders), 2);

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    amp_name = [subj_name,'_amp'];
    
    figure;
    
    h_all = nan(50, 2*no_periods, 2);
    
    b_all = nan(50, 2*no_periods, 2);
    
    for pd = 1:no_periods
        
        A_pd = A(periods(pd, 1):periods(pd, 2), :, 3);
        
        A_tot_pd = sum(A(periods(pd, 1):periods(pd, 2), :, :), 3);
        
        A_pd = A_pd./A_tot_pd;
        
        A_rk = tiedrank(A_pd);
        
        fprintf(fid_all, '%f\t%f\t%f\t%f\t%f\n', [pd*ones(size(A_pd), 1) A_pd A_rk*diag(1./max(A_rk))]');
        
        %% Plotting amplitude of channel 1 by amplitude of channel 2.
        
        subplot(3, no_periods, pd)
        
        r = nancorr(A_rk(:, 1), A_rk(:, 2));
        
        [histogram, centers] = hist3(A_rk, [50 50]);
        
        colormap('default')
        
        imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
        
        h = colorbar;
        
        title({[folder, ' ', period_label{pd}]; ['Correlation = ', num2str(r)]})
        
        ylabel(h, 'Prop. Observed')
        
        axis xy
        
        xlabel([chan_labels{1},' Power (Rank)'])
        
        ylabel([chan_labels{2}, ' Power (Rank)'])
        
        %% Computing distribution of power in one channel for quartiles of power in other channel.
        
        for ch1 = 1:2
            
            top_q_cutoff = quantile(A_pd(:, ch1), 0.75);
            
            top_q_indicator = A_pd(:, ch1) >= top_q_cutoff;
            
            other_amp = A_pd(top_q_indicator, 3 - ch1);
            
            [h, b] = hist(other_amp, 50);
            
            h_all(:, pd, ch1) = h/sum(h);
            
            b_all(:, pd, ch1) = b;
            
            q_legend{pd} = [period_label{pd}, ', Top 25%']; 
            
            bot_q_cutoff = quantile(A_pd(:, ch1), 0.25);
            
            bot_q_indicator = A_pd(:, ch1) <= bot_q_cutoff;
            
            other_amp = A_pd(bot_q_indicator, 3 - ch1);
            
            [h, b] = hist(other_amp, 50);
            
            h_all(:, 2 + pd, ch1) = h/sum(h);
            
            b_all(:, 2 + pd, ch1) = b;
            
            q_legend{2 + pd} = [period_label{pd}, ', Bottom 25%']; 
            
        end
        
        %% Computing covariance of power in two channels.
        
        [A_corr, lags] = xcorr(detrend(A_pd(:, 1), 'constant'), detrend(A_pd(:, 2), 'constant'));
        
        % A_corr = A_corr/prod(std(A_pd));
        
        subplot(3, no_periods, 2*no_periods + pd)
        
        lag_indices = abs(lags) <= win_size;
        
        lags = lags(lag_indices);
        
        A_corr = A_corr(lag_indices);
        
        plot(lags, A_corr);
        
        hold on
        
        all_corr(fo, :, pd) = A_corr;
        
        [~, max_index] = max(A_corr);
        
        max_lag = lags(max_index);
        
        all_lags(fo, pd) = max_lag;
        
        plot([max_lag; max_lag], [min(A_corr); max(A_corr)], 'r')
        
        axis tight
        
        title([folder,', Cross-Correlation of Amplitude, ', period_label{pd}])
        
        xlabel('Lag (ms)')
        
    end
    
    %% Plotting distributions computed above.
    
    for ch1 = 1:2
        
        subplot(3, no_periods, no_periods + ch1)
        
        plot(b_all(:, :, ch1), h_all(:, :, ch1))
        
        legend(q_legend)
        
        title([chan_labels{3 - ch1}, ' Amp. Dist. for Top/Bottom Quartile of ', chan_labels{ch1}, ' Power'])
        
        ylabel('Proportion of Observations')
        
        xlabel([chan_labels{3 - ch1}, ' Beta Amp.'])
        
    end
    
    save_as_pdf(gcf, amp_name)
    
end

fclose(fid_all);

% end

all_amp_data = load(all_amp_name);

all_pds = all_amp_data(:, 1);

norm_labels = {'Power','Power (Rank)'};

no_norms = length(norm_labels);

figure

for norm = 1:no_norms
    
    for pd = 1:no_periods
        
        %% Plotting amplitude against amplitude.
        
        subplot(no_norms + 2, no_periods, (norm - 1)*no_periods + pd)
        
        A_pd = all_amp_data(all_pds == pd, (norm - 1)*no_periods + (2:3));
        
        r = nancorr(A_pd(:, 1), A_pd(:, 2));
        
        [histogram, centers] = hist3(A_pd, [50 50]);
        
        colormap('default')
        
        imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
        
        h = colorbar;
        
        title({[chan_labels{1},' by ',chan_labels{2},', ',period_label{pd}];['Correlation = ',num2str(r)]})
        
        ylabel(h, 'Prop. Observed')
        
        axis xy
        
        xlabel([chan_labels{1}, ' ', norm_labels{norm}])
        
        ylabel([chan_labels{2}, ' ', norm_labels{norm}])
        
        if norm == 1
            
            %% Computing distribution of ... by quartile of ...
        
            for ch1 = 1:2
                
                top_q_cutoff = quantile(A_pd(:, ch1), 0.75);
                
                top_q_indicator = A_pd(:, ch1) >= top_q_cutoff;
                
                other_amp = A_pd(top_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, pd, ch1) = h/sum(h);
                
                b_all(:, pd, ch1) = b;
                
                q_legend{pd} = [period_label{pd}, ', Top 25%'];
                
                bot_q_cutoff = quantile(A_pd(:, ch1), 0.25);
                
                bot_q_indicator = A_pd(:, ch1) <= bot_q_cutoff;
                
                other_amp = A_pd(bot_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, 2 + pd, ch1) = h/sum(h);
                
                b_all(:, 2 + pd, ch1) = b;
                
                q_legend{2 + pd} = [period_label{pd}, ', Bottom 25%'];
                
            end
            
        end
        
    end
    
    if norm == 1
        
        %% Plotting distribution of ... by quartile of ...
        
        for ch1 = 1:2
            
            subplot(no_norms + 2, no_periods, no_norms*no_periods + ch1)
            
            plot(b_all(:, :, ch1), h_all(:, :, ch1))
            
            legend(q_legend)
            
            title([chan_labels{3 - ch1}, ' Amp. Dist. for Top/Bottom Quartile of ', chan_labels{ch1}, ' Amp.'])
            
            ylabel('Proportion of Observations')
            
            xlabel([chan_labels{3 - ch1}, ' Beta Amp.'])
            
        end
        
        %% Plotting Cross-Correlation.
        
        subplot(no_norms + 2, no_periods, (no_norms + 2)*no_periods - 1)
        
        plot(-win_size:win_size, reshape(mean(all_corr), length(-win_size:win_size), 2))
        
        axis tight
        
        title(['Mean Cross-Corr., ', chan_labels{1}, ' Amp. by ', chan_labels{2}, ' Amp.'])
        
        xlabel('Lag (ms)')
        
        ylabel('Cross-Corr.')
        
        subplot(no_norms + 2, no_periods, (no_norms + 2)*no_periods)
        
        plot([1 2], all_lags', '-*')
        
        set(gca, 'XTick', [1 2], 'XTickLabel', period_label)
        
        xlim([.75 2.25])
        
        ylim([-15, 15])
        
        title('Time Lag, Pre- and Post-Infusion')
        
        ylabel('Time Lag (ms)')
        
    end
    
end

save_as_pdf(gcf, all_amp_name(1:end-4))