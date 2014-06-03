function PD_beta_epochs(sd_lim, win_size, step_size)

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

% basetimes = [300 1200 1800 600 1800];
% 
% infusetimes = [390 240 300 450 510];

% bands = [1 4; 4 12; 15 30; 30 60; 60 90; 90 110; 120 180];
% band_names = {'delta','theta','beta','lgamma','hgamma','HFO'};
% no_bands = length(band_names);

sampling_freq = 1000;

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    subj_name = [folder,'/',prefix];
    
    load([subj_name,'_all_channel_data_dec.mat']) 
    
    load([subj_name,'_all_channel_data_dec_HAP.mat'])
    
    t = (1:size(A,1))/sampling_freq;
    
    beta_amp = A(:,:,3);
    
%     beta_cutoff = mean(beta_amp) + sd_lim*std(beta_amp);
    
    figure;
    
    for ch = 1:2
         
        beta_listname = [subj_name,'_ch',num2str(ch), '_beta.list'];
        
        fid_list = fopen(beta_listname, 'w');
        
        fid_P_list = fopen([beta_listname(1:end-5), 'P.list'],'w');
        
        fid_win_list = fopen([beta_listname(1:end-5),'_win.list'],'w');
       
        ba_flipped = [flipud(beta_amp(1:10*sampling_freq,ch)); beta_amp(:,ch); flipud(beta_amp((end-10*sampling_freq+1):end,ch))];
        
        ba_conv = conv(ba_flipped,ones(20*sampling_freq,1)/(20*sampling_freq),'same');
        
        ba_smooth = ba_conv((10*sampling_freq+1):(end-10*sampling_freq));
        
%         beta_high = beta_amp(:,ch) >= beta_cutoff(ch);
        
        beta_cutoff = mean(ba_smooth) + sd_lim*std(ba_smooth);
        
        beta_high = ba_smooth >= beta_cutoff;
       
        dbh = diff(beta_high);
        
        beta_start = find(dbh == 1) + 1;
        
        beta_end = find(dbh == -1);
        
        if beta_end(1) < beta_start(1)
            
            beta_start = [1; beta_start];
        
        end
            
        if beta_start(end) > beta_end(end)
            
            beta_end = [beta_end; size(beta_amp,1)];
            
        end
        
        beta_blocks = [beta_start beta_end];
        
        beta_lengths = diff(beta_blocks,[],2) + 1;
        
        beta_blocks(beta_lengths < win_size, :) = [];

        beta_lengths = diff(beta_blocks,[],2) + 1;
        
        subplot(2,1,ch)
        
%         plot(t,beta_amp(:,ch),'k')
        
        plot(t,beta_amp(:,ch),'k',t,ba_smooth,'b')

        hold on
        
        for b = 1:size(beta_blocks,1)
            
            plot(t(beta_blocks(b,1):beta_blocks(b,2)),beta_amp(beta_blocks(b,1):beta_blocks(b,2),ch),'g')
            
            plot(t(beta_blocks(b,1):beta_blocks(b,2)),ba_smooth(beta_blocks(b,1):beta_blocks(b,2)),'r')
            
            beta_name = [subj_name,'_ch',num2str(ch),'_beta_block',num2str(b)];
            
            no_epochs = floor((beta_lengths(b) - win_size)/step_size) + 1;
        
            for e = 1:no_epochs
               
                epoch_name = [beta_name,'_epoch',num2str(e),'.txt'];
                
                P_name = [beta_name,'_epoch',num2str(e),'_P.txt'];
                
                epoch_start = beta_blocks(b,1) + e*step_size;
                
                epoch_end = beta_blocks(b,1) + e*step_size + win_size;
                
                fid = fopen(epoch_name,'w');
                
                fprintf(fid, '%f\t%f\n', PD_dec(epoch_start:epoch_end, :)');
                
                fclose(fid);
                
                fid = fopen(P_name, 'w');
                
                fprintf(fid, '%f\t%f\n', P(epoch_start:epoch_end, :, 3)');
                
                fclose(fid);
                
                fprintf(fid_list, '%s\n', epoch_name);
                
                fprintf(fid_P_list, '%s\n', P_name);
                
                fprintf(fid_win_list, '%d\n', round(mean([epoch_start; epoch_end])));
                
            end
            
        end
        
    end
    
    saveas(gcf,[subj_name,'_beta_',num2str(sd_lim),'sd_',num2str(win_size),'win_',num2str(step_size),'step.fig'])
          
end