function PD_beta_epochs_rel_infusion(subject_mat, outlier_lim, sd_lim, win_size)

load(subject_mat)

sampling_freq = 1000;

pd_label = {'pre', 'post'};

period_label = {'Pre-Infusion','Post-Infusion'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
        
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    t = (1:size(A,1))/sampling_freq;
    
    beta_amp = A(:,:,3); ba_smooth = nan(size(beta_amp)); beta_high = nan(size(beta_amp));
    
    pd_limits = [1 base_index; (base_index + 1) min(length(t), base_index + 1500*sampling_freq)];
    
    % beta_cutoff = mean(beta_amp(1:base_index, :)) + sd_lim*std(beta_amp(1:base_index, :));
    
    figure;
    
    for ch = 1:2
        
        norm_data = (PD_dec(:, ch) - mean(PD_dec(:, ch)))/std(PD_dec(:, ch));
        
        ba_flipped = [flipud(beta_amp(1:win_size, ch)); beta_amp(:,ch); flipud(beta_amp((end-win_size+1):end,ch))];
        
        ba_conv = conv(ba_flipped, ones(win_size, 1)/win_size, 'same');
        
        ba_smooth(:, ch) = ba_conv((win_size+1):(end-win_size));
        
        beta_cutoff = mean(ba_smooth(1:base_index, ch)) + sd_lim*std(ba_smooth(1:base_index, ch));
        
        beta_high(:, ch) = ba_smooth(:, ch) >= beta_cutoff;
        
        for pd = 1:size(pd_limits,1)
            
            beta_listname = [subj_name,'_ch',num2str(ch),'_beta_',pd_label{pd},'.list'];
            
            fid_list = fopen(beta_listname, 'w');
            
            fid_P_list = fopen([beta_listname(1:end-5), '_P.list'],'w');
            
            fid_win_list = fopen([beta_listname(1:end-5), '_win.list'],'w');
            
            subplot(2, 2, (ch-1)*2 + pd)
            
            % plot(t,beta_amp(:,ch),'k')
            
            title([chan_labels{ch}, ' Beta Segments ', period_label{pd}])
            
            plot(t(pd_limits(pd,1):pd_limits(pd,2)), beta_amp(pd_limits(pd,1):pd_limits(pd,2), ch), 'k',...
                t(pd_limits(pd,1):pd_limits(pd,2)), ba_smooth(pd_limits(pd,1):pd_limits(pd,2), ch), 'b')
            
            axis tight
            
            hold on
            
            % ba_pd = ba_smooth(pd_limits(pd,1):pd_limits(pd,2));
            
            % beta_cutoff = mean(ba_pd) + sd_lim*std(ba_pd);
            
            % beta_high = ba_pd >= beta_cutoff;
            
            bh_pd = beta_high(pd_limits(pd,1):pd_limits(pd,2), ch);
            
            dbh = diff(bh_pd);%beta_high);
            
            beta_start = find(dbh == 1) + 1 + pd_limits(pd,1) - 1;
            
            beta_end = find(dbh == -1) + pd_limits(pd,1) - 1;
            
            if ~isempty(beta_end) && ~isempty(beta_start)
                
                if beta_end(1) < beta_start(1)
                    
                    beta_start = [pd_limits(pd,1); beta_start];
                    
                end
                
                if beta_start(end) > beta_end(end)
                    
                    beta_end = [beta_end; pd_limits(pd,2)];
                    
                end
                
                beta_blocks = [beta_start beta_end];
                
                beta_lengths = diff(beta_blocks,[],2) + 1;
                
                beta_blocks(beta_lengths < win_size, :) = [];
                
                beta_out = [];
                
                index = 1;
                
                for b = 1:size(beta_blocks,1)
                        
                    beta_name = [subj_name,'_ch',num2str(ch),'_beta_',pd_label{pd},'_block',num2str(b),'.txt'];
                    
                    if any(abs(norm_data(beta_blocks(b,1):beta_blocks(b,2))) > outlier_lim)
                        
                        beta_out(index) = b;
                        
                        plot(t(beta_blocks(b,1):beta_blocks(b,2)), beta_amp(beta_blocks(b,1):beta_blocks(b,2), ch), 'c')
                        
                        plot(t(beta_blocks(b,1):beta_blocks(b,2)), ba_smooth(beta_blocks(b,1):beta_blocks(b,2), ch), 'm')
                        
                        index = index + 1;
                        
                    else
                        
                        plot(t(beta_blocks(b,1):beta_blocks(b,2)), beta_amp(beta_blocks(b,1):beta_blocks(b,2), ch), 'g')
                        
                        plot(t(beta_blocks(b,1):beta_blocks(b,2)), ba_smooth(beta_blocks(b,1):beta_blocks(b,2), ch), 'r')
                        
                        P_name = [beta_name(1:end-4),'_P.txt'];
                        
                        fid = fopen(beta_name,'w');
                        
                        fprintf(fid, '%f\t%f\n', PD_dec(beta_blocks(b,1):beta_blocks(b,2), :)');
                        
                        fclose(fid);
                        
                        fid = fopen(P_name, 'w');
                        
                        fprintf(fid, '%f\t%f\n', P(beta_blocks(b,1):beta_blocks(b,2), :, 3)');
                        
                        fclose(fid);
                        
                        fprintf(fid_list, '%s\n', beta_name);
                        
                        fprintf(fid_P_list, '%s\n', P_name);
                        
                        fprintf(fid_win_list, '%d\t%d\t%d\t%f\t%f\n', b, beta_blocks(b,:), median(beta_amp(beta_start:beta_end, :)));
                        
                    end
                    
                end
                    
                no_outs = length(beta_out);
                
                if no_outs > 0
                    
                    figure;
                    
                    [s_r, s_c] = subplot_size(no_outs);
                    
                    for o = 1:no_outs
                        
                        bo = beta_out(o);
                        
                        subplot(s_r, s_c, o)
                        
                        plot(t(beta_blocks(bo,1):beta_blocks(bo,2)), norm_data(beta_blocks(bo,1):beta_blocks(bo,2)))
                        
                        axis tight
                        
                    end
                    
                    save_as_pdf(gcf, [beta_listname(1:end-5),'_out'])
                    
                    close(gcf)
                    
                end
                
            end
            
        end
        
    end
    
    save_as_pdf(gcf,[subj_name,'_beta_',pd_label{pd},'_',num2str(outlier_lim),'out_',num2str(sd_lim),'sd_',num2str(win_size),'win'])
          
end