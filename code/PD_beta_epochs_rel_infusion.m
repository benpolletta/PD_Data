function PD_beta_epochs_rel_infusion(sd_lim, win_size)

load('subjects.mat')

sampling_freq = 1000;

pd_label = {'pre', 'post'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    base_index = basetimes(fo)*sampling_freq;
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat'])
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    t = (1:size(A,1))/sampling_freq;
    
    beta_amp = A(:,:,3);
    
    %     beta_cutoff = mean(beta_amp) + sd_lim*std(beta_amp);
    
    pd_limits = [1 base_index; (base_index + 1) length(t)];
    
    figure;
    
    for ch = 1:2
        
        ba_flipped = [flipud(beta_amp(1:10*sampling_freq,ch)); beta_amp(:,ch); flipud(beta_amp((end-10*sampling_freq+1):end,ch))];
        
        ba_conv = conv(ba_flipped,ones(20*sampling_freq,1)/(20*sampling_freq),'same');
        
        ba_smooth = ba_conv((10*sampling_freq+1):(end-10*sampling_freq));
        
        %         beta_high = beta_amp(:,ch) >= beta_cutoff(ch);
        
        for pd = 1:size(pd_limits,1)
            
            beta_listname = [subj_name,'_ch',num2str(ch),'_beta_',pd_label{pd},'.list'];
            
            fid_list = fopen(beta_listname, 'w');
            
            fid_P_list = fopen([beta_listname(1:end-5), '_P.list'],'w');
            
            fid_win_list = fopen([beta_listname(1:end-5),'_win.list'],'w');
            
            subplot(2,2,(ch-1)*2 + pd)
            
            %         plot(t,beta_amp(:,ch),'k')
            
            plot(t,beta_amp(:,ch),'k',t,ba_smooth,'b')
            
            hold on
            
            ba_pd = ba_smooth(pd_limits(pd,1):pd_limits(pd,2));
            
            beta_cutoff = mean(ba_pd) + sd_lim*std(ba_pd);
            
            beta_high = ba_pd >= beta_cutoff;
            
            dbh = diff(beta_high);
            
            beta_start = find(dbh == 1) + 1 + pd_limits(pd,1) - 1;
            
            beta_end = find(dbh == -1) + pd_limits(pd,1) - 1;
            
            if beta_end(1) < beta_start(1)
                
                beta_start = [pd_limits(pd,1); beta_start];
                
            end
            
            if beta_start(end) > beta_end(end)
                
                beta_end = [beta_end; pd_limits(pd,2)];
                
            end
            
            beta_blocks = [beta_start beta_end];
            
            beta_lengths = diff(beta_blocks,[],2) + 1;
            
            beta_blocks(beta_lengths < win_size, :) = [];
            
            for b = 1:size(beta_blocks,1)
                
                plot(t(beta_blocks(b,1):beta_blocks(b,2)),beta_amp(beta_blocks(b,1):beta_blocks(b,2),ch),'g')
                
                plot(t(beta_blocks(b,1):beta_blocks(b,2)),ba_smooth(beta_blocks(b,1):beta_blocks(b,2)),'r')
                
                beta_name = [subj_name,'_ch',num2str(ch),'_beta_',pd_label{pd},'_block',num2str(b),'.txt'];
                
                P_name = [beta_name(1:end-4),'_P.txt'];
                
                fid = fopen(beta_name,'w');
                
                fprintf(fid, '%f\t%f\n', PD_dec(beta_blocks(b,1):beta_blocks(b,2), :)');
                
                fclose(fid);
                
                fid = fopen(P_name, 'w');
                
                fprintf(fid, '%f\t%f\n', P(beta_blocks(b,1):beta_blocks(b,2), :, 3)');
                
                fclose(fid);
                
                fprintf(fid_list, '%s\n', beta_name);
                
                fprintf(fid_P_list, '%s\n', P_name);
                
                fprintf(fid_win_list, '%d\t%d\t%d\t%f\t%f\n', b, beta_blocks(b,:), median(beta_amp(beta_start:beta_end,:)));
                
            end
            
        end
        
    end
    
    saveas(gcf,[subj_name,'_beta_',pd_label{pd},'_',num2str(sd_lim),'sd_',num2str(win_size),'win.fig'])
          
end