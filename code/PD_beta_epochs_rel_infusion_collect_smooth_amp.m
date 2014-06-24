function PD_beta_epochs_rel_infusion_collect_smooth_amp(subject_mat, outlier_lim, sd_lim, win_size, smooth_size)

par_name = [num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(smooth_size),'smooth'];

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
        
        h_all = nan(50, no_periods, 2);
        
        b_all = nan(50, no_periods, 2);
        
        for pd = 1:no_periods
            
            beta_listname = [beta_name,'_',pd_label{pd},'_',par_name,'.list'];
            
            % Loading block numbers, epoch start and end indices.
            [blocks, beta_starts, beta_ends] = text_read([beta_listname(1:end-5),'_win.list'],'%f%f%f%*[^\n]');
            
            beta_amp_name = [beta_listname(1:end-5),'_amp.txt'];
            
            % if isempty(dir(beta_amp_name))
                
                fid = fopen(beta_amp_name, 'w');
                
                no_blocks = length(blocks);
                
                %% Computing smoothed frequency and phase difference for each block.
                
                for b = 1:no_blocks
                    
                    A_block = A(beta_starts:beta_ends, :, 3);
                    
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
                    
                end
                
                fclose(fid);
                
            % end
            
            all_beta_data = load(beta_amp_name);
                
            fprintf(fid_mat(ch), '%f\t%f\t%f\t%f\n', [pd*ones(size(all_beta_data, 1), 1) all_beta_data]');
            
            if ~isempty(all_beta_data)
                
                all_A = all_beta_data(:, 2:3);
                
            else
                
                all_A = [nan nan];
                
            end
            
            subplot(2, 2, pd)
            
            r = nancorr(all_A(:, 1), all_A(:, 2));
            
            [histogram, centers] = hist3(all_A, [50 50]);
            
            colormap('default')
            
            imagesc(centers{1}, centers{2}, histogram/sum(sum(histogram)))
            
            h = colorbar;
            
            title({[folder,' ',chan_labels{ch},' High Beta ',period_label{pd}];['Correlation = ',num2str(r)]})
            
            ylabel(h, 'Prop. Observed')
            
            axis xy
            
            xlabel([chan_labels{1},' Power'])
            
            ylabel([chan_labels{2}, ' Power'])
            
            for ch1 = 1:2
            
                top_q_cutoff = quantile(all_A(:, ch1), 0.75);
                
                top_q_indicator = all_A(:, ch1) >= top_q_cutoff;
                
                other_amp = all_A(top_q_indicator, 3 - ch1);
                
                [h, b] = hist(other_amp, 50);
                
                h_all(:, pd, ch1) = h/sum(h);
                
                b_all(:, pd, ch1) = b;
                
            end
            
        end
        
        for ch1 = 1:2
            
            subplot(2, 2, 2 + ch1)
            
            plot(b_all(:, :, ch1), h_all(:, :, ch1))
            
            legend(period_label)
            
            title([chan_labels{3 - ch1}, ' Amp. Dist. for Top Quartile of ', chan_labels{ch1}, ' Power'])
            
            ylabel('Proportion of Observations')
            
            xlabel([chan_labels{3 - ch1}, ' Beta Amp.'])
            
        end
        
        save_as_pdf(gcf,[beta_name,'_ri_aa'])
        
    end
    
end

close('all')

for ch = 1:no_channels
    
    fclose(fid_mat(ch));
    
end