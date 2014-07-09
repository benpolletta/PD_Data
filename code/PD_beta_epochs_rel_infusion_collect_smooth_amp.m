function PD_beta_epochs_rel_infusion_collect_smooth_amp(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

colors = {'b', 'g', 'r'};

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

cov_length = 2*win_size + 1;

cov_format = make_format(cov_length, 'f');

load(subject_mat)

sampling_freq = 1000;

pd_label = {'pre','post'};

period_label = {'Pre-Infusion','Post-Infusion'};

chan_labels = {chan_labels{:}, 'Both', [chan_labels{1},' Not ',chan_labels{2}], [chan_labels{2},' Not ',chan_labels{1}], [chan_labels{1},' High ',chan_labels{2},' Low '], [chan_labels{1},' Low ',chan_labels{2},' High']};

ch_label = {'ch1', 'ch2', 'ch1_ch2', 'ch1_nch2', 'ch2_nch1', 'ch1_lch2', 'ch2_lch1'};

no_channels = length(ch_label);

fid_mat = nan(no_channels,1);

all_beta_name = cell(no_channels,1);

for ch = 1:no_channels
    
    all_beta_name{ch} = [subject_mat(1:(end-length('_subjects.mat'))),'_',par_name,'_beta_',ch_label{ch}];
    
    fid_mat(ch) = fopen([all_beta_name{ch},'_amp.txt'], 'w');
    
    cov_fid_mat(ch) = fopen([all_beta_name{ch},'_cov.txt'], 'w');
    
    lag_fid_mat(ch) = fopen([all_beta_name{ch},'_lag.txt'], 'w');
    
end

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};  
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    periods = [1 base_index; (base_index + 1) size(A,1)];
    
    no_periods = size(periods,1);
    
    for ch = 1:no_channels
        
        beta_name = [subj_name,'_',ch_label{ch},'_beta'];
        
        figure;
        
        % h_all = nan(50, no_periods, 2);
        
        % b_all = nan(50, no_periods, 2);
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            beta_amp_name = [beta_listname(1:end-5),'_amp.txt'];
            
            beta_cov_name = [beta_listname(1:end-5),'_cov.txt'];
            
            beta_lag_name = [beta_listname(1:end-5),'_lag.txt'];
            
            % if isempty(dir(beta_amp_name))
                
                fid = fopen(beta_amp_name, 'w');
                
                cov_fid = fopen(beta_cov_name, 'w');
                
                lag_fid = fopen(beta_lag_name, 'w');
                
                no_blocks = length(blocks);
                
                %% Computing amplitude and covariance sequence for each block.
                
                for b = 1:no_blocks
                    
                    A_block = A(beta_starts(b):beta_ends(b), :, 3);
                    
                    % A_smooth = nan(size(A_block, 1), 2);
                    % 
                    % for ch1 = 1:2
                    % 
                    %     A_flipped = [flipud(A_block(1:smooth_winsize, ch1)); A_block(:, ch1); flipud(A_block((end-smooth_winsize+1):end, ch1))];
                    % 
                    %     A_conv = conv(A_flipped, hann(smooth_winsize)/sum(hann(smooth_winsize)), 'same');
                    %     % A_conv = conv(A_flipped, ones(smooth_winsize,1)/smooth_winsize, 'same');
                    % 
                    %     A_smooth(:, ch1) = A_conv(smooth_winsize + (1:size(A_block, 1)));
                    % 
                    % end
                    
                    fprintf(fid, '%f\t%f\t%f\n', [b*ones(size(A_block, 1), 1) A_block]');
                    
                    % Computing and printing covariance.
                    
                    [A_cov, lags] = xcorr(detrend(A_block(:, 1), 'constant'), detrend(A_block(:, 2), 'constant'));
                    
                    % A_cov = A_cov/prod(var(A_block));
                    
                    A_cov = A_cov(abs(lags) <= win_size);
                    
                    lags = lags(abs(lags) <= win_size);
                    
                    fprintf(cov_fid, cov_format, A_cov);
                    
                    % Computing and printing lag with max correlation.
                    
                    [max_cov, max_index] = max(A_cov);
                    
                    max_lag = lags(max_index);
                    
                    fprintf(lag_fid, '%f\t%f\n', max_lag, max_cov);
                    
                end
                
                fclose(fid);
                
            % end
           
            %% Printing amplitude for all beta blocks, for each channel.
            
            all_beta_data = load(beta_amp_name);
                
            fprintf(fid_mat(ch), '%f\t%f\t%f\t%f\n', [pd*ones(size(all_beta_data, 1), 1) all_beta_data]');
            
            all_cov_data = load(beta_cov_name);
                
            fprintf(cov_fid_mat(ch), ['%f\t', cov_format], [pd*ones(size(all_cov_data, 1), 1) all_cov_data]');
            
            all_lag_data = load(beta_lag_name);
                
            fprintf(lag_fid_mat(ch), '%f\t%f\t%f\n', [pd*ones(size(all_lag_data, 1), 1) all_lag_data]');
            
            %% Plotting amplitude correlation for all beta blocks, for each channel.
            
            if ~isempty(all_cov_data)
            
                all_cov = all_cov_data;
            
            else
            
                all_cov = nan(1, cov_length);
            
            end
            
            subplot(1, 2, pd)
            
            plot(-win_size:win_size, mean(all_cov)')
            
            axis tight
            
            box off
            
            % %% Plotting amplitude correlation for all beta blocks, for each channel.
            % 
            % if ~isempty(all_beta_data)
            % 
            %     all_A = all_beta_data(:, 2:3);
            % 
            % else
            % 
            %     all_A = [nan nan];
            % 
            % end
            % 
            % subplot(2, 2, pd)
            % 
            % r = nancorr(all_A(:, 1), all_A(:, 2));
            % 
            % [histogram, centers] = hist3(all_A, [50 50]);
            % 
            % colormap('default')
            % 
            % imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
            % 
            % h = colorbar;
            % 
            % title({[folder,' ',chan_labels{ch},' High Beta ',period_label{pd}];['Correlation = ',num2str(r)]})
            % 
            % ylabel(h, 'Prop. Observed')
            % 
            % axis xy
            % 
            % xlabel([chan_labels{1},' Power'])
            % 
            % ylabel([chan_labels{2}, ' Power'])
            
            % %% Calculating distibution of channel 2 by quartile of channel 1, and vice-versa.
            %
            % for ch1 = 1:2
            % 
            %     top_q_cutoff = quantile(all_A(:, ch1), 0.75);
            % 
            %     top_q_indicator = all_A(:, ch1) >= top_q_cutoff;
            % 
            %     other_amp = all_A(top_q_indicator, 3 - ch1);
            % 
            %     [h, b] = hist(other_amp, 50);
            % 
            %     h_all(:, pd, ch1) = h/sum(h);
            % 
            %     b_all(:, pd, ch1) = b;
            % 
            % end
            
        end
        
        % %% Plotting distribution of amplitude for top and bottom quartile in other channel.
        %
        % for ch1 = 1:2
        % 
        %     subplot(2, 2, 2 + ch1)
        % 
        %     plot(b_all(:, :, ch1), h_all(:, :, ch1))
        % 
        %     legend(period_label)
        % 
        %     title([chan_labels{3 - ch1}, ' Amp. Dist. for Top Quartile of ', chan_labels{ch1}, ' Power'])
        % 
        %     ylabel('Proportion of Observations')
        % 
        %     xlabel([chan_labels{3 - ch1}, ' Beta Amp.'])
        % 
        % end
        % 
        % save_as_pdf(gcf, [beta_name, '_ri_aa'])
        
        save_as_pdf(gcf, [beta_name, '_ri_a_cov'])
        
        close(gcf)
        
    end
    
end

close('all')

for ch = 1:no_channels
    
    fclose(fid_mat(ch)); fclose(cov_fid_mat(ch)); fclose(lag_fid_mat(ch));
    
    figure
    
    chan_cov = load([all_beta_name{ch}, '_cov.txt']);
    
    if ~isempty(chan_cov)
        
        chan_pd = chan_cov(:, 1);
        
        chan_cov = chan_cov(:, 2:end);
        
        mean_cov = nan(size(chan_cov, 2), 2); std_cov = nan(size(chan_cov, 2), 2);
        
        for pd = 1:2
            
            mean_cov(:, pd) = mean(chan_cov(chan_pd == pd, :))';
            
            std_cov(:, pd) = std(chan_cov(chan_pd == pd, :))';
            
        end
        
        subplot(1, 3, 1)
        
        h = plot(-win_size:win_size, mean_cov);
        
        hold on
        
        plot(-win_size:win_size, mean_cov + std_cov, '--')
        
        plot(-win_size:win_size, mean_cov - std_cov, '--')
        
        axis tight
        
        title({chan_labels{ch};'Mean \pm S.D. Cross-Correlation'})
        
        legend(h, period_label)
        
    end
    
    chan_lag = load([all_beta_name{ch}, '_lag.txt']);
        
    if ~isempty(chan_lag)
        
        chan_pds = chan_lag(:, 1);
        
        chan_lag = chan_lag(:, 2:end);
        
        for pd = 1:2
        
            subplot(1, 3, 2)
            
            scatter(chan_lag(chan_pds == pd, 1), chan_lag(chan_pds == pd, 2), colors{pd})
        
            hold on
    
            [h, b] = hist(chan_lag(chan_pds == pd, 1), -100:10:100);
            
            histograms(pd, :) = h/sum(h);
            
            bins(pd, :) = b;
            
            axis tight
            
            title({chan_labels{ch};'Max. Cross-Corr. by Lag Time of Max. Cross-Corr.'})
            
        end
    
        subplot(1, 3, 3)
        
        plot(bins', histograms')
        
        title({chan_labels{ch};'Histogram of Lag for Max. Cross-Corr.'})
        
        xlabel('Lag (ms)')
        
    end
    
    save_as_pdf(gcf, [all_beta_name{ch}, '_cov'])
    
end