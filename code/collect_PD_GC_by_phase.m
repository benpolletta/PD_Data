function collect_PD_GC_by_phase

phase_bin_edges = -pi:pi/12:pi;

phase_centers = -mean([phase_bin_edges(1:end-1); phase_bin_edges(2:end)])';

c_order = [1 0 1; 0 1 1];

chan_legend = {'Striatal --> Motor','Motor --> Striatal'};

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    cd (folder)
    
    master_listname = [folder,'_all_channel_data_dec_phase_master.list'];
    
    master_list = text_read(master_listname,'%s%*[^\n]');
    
    no_lists = length(master_list);
   
    phases = [];
    
    GC_collected = [];
    
    p_collected = [];
    
    GC_index = 0;
   
    mean_GC = nan(no_lists,2);
    
    std_GC = nan(no_lists,2);
    
    mean_p = nan(no_lists,2);
    
    std_p = nan(no_lists,2);
    
    sig_count = nan(no_lists,2);
    
    for l = 1:no_lists
        
        listname = master_list{l};
        
        load([listname(1:end-5),'_GC.mat']) % Loading 'Orders', 'Errors', 'All_GC', 'All_pval', 'All_sig'.
        
        no_epochs = size(All_GC,1);
        
        if no_epochs >= 1
            
            phases(GC_index+(1:no_epochs),1) = phase_centers(l);
            
            GC_collected(GC_index+(1:no_epochs),:) = All_GC;
            
            p_collected(GC_index+(1:no_epochs),:) = All_pval;
            
        end
        
        GC_index = GC_index + no_epochs;
        
        mean_GC(l,:) = nanmean(All_GC);
        
        std_GC(l,:) = nanstd(All_GC);
        
        mean_p(l,:) = nanmean(1-All_pval);
        
        std_p(l,:) = nanstd(1-All_pval);
        
        sig_count(l,:) = sum(All_sig == 1);
        
    end
    
    save([master_listname(1:end-5),'_GC.mat'],'phases','GC_collected','mean_GC','std_GC','mean_p','std_p','sig_count')
    
    mean_vec = sum(mean_GC.*repmat(exp(sqrt(-1)*phase_centers),1,2));
    
    p_vec = sum(mean_p.*repmat(exp(sqrt(-1)*phase_centers),1,2));
    
    sig_vec = sum(sig_count.*repmat(exp(sqrt(-1)*phase_centers),1,2));
    
    phase_for_plot = [phase_centers; phase_centers(1)];
    
    mean_GC = [mean_GC; mean_GC(1,:)];
    
    std_GC = [std_GC; std_GC(1,:)];
    
    mean_p = [mean_p; mean_p(1,:)];
    
    std_p = [std_p; std_p(1,:)];
    
    sig_count = [sig_count; sig_count(1,:)];
    
    figure; % These plots could be folded into a for loop at some point, if you have time on your hands.
    
    subplot(1,2,1)
    
    polar_lim=polar(0, max(max(max(mean_GC + std_GC)), max(abs(mean_vec))), '.'); % Dummy plot to set radial axis length.

    set(polar_lim,'Marker','none'); % Dummy plot disappears.
    
    hold on
    
    for ch = 1:2
        
        h1(ch) = polar(phase_for_plot, mean_GC(:,ch), '-');
        
        set(h1(ch),'Color',c_order(ch,:),'LineWidth',2)
        
        h = polar(phase_for_plot, mean_GC(:,ch) + std_GC(:,ch), ':');
        
        set(h,'Color',c_order(ch,:),'LineWidth',2)
        
        h = compass(mean_vec(ch));
        
        set(h,'Color',c_order(ch,:),'LineWidth',2)
        
    end
    
    title([folder,', GC (Mean \pm SD) by Phase Diff. (S - M)'])
    
%     colormap(c_order)
%     
%     colorbar('YTick',(1:no_f_bins)+.5,'YTickLabel',freq_legend)
    
    legend(h1,chan_legend)
    
    subplot(1,2,2)
    
    polar_lim=polar(0, max(max(max(mean_p + std_p)), max(abs(mean_vec))), '.'); % Dummy plot to set radial axis length.

    set(polar_lim,'Marker','none'); % Dummy plot disappears.
    
    hold on
    
    for ch = 1:2
        
        h1(ch) = polar(phase_for_plot, mean_p(:,ch), '-');
        
        set(h1(ch),'Color',c_order(ch,:),'LineWidth',2)
        
        h = polar(phase_for_plot, mean_p(:,ch) + std_p(:,ch), ':');
        
        set(h,'Color',c_order(ch,:),'LineWidth',2)
        
        h = compass(mean_vec(ch));
        
        set(h,'Color',c_order(ch,:),'LineWidth',2)
        
    end
    
    title([folder,', GC p-Value (Mean \pm SD) by Phase Diff. (S - M)'])
    
%     colormap(c_order)
%     
%     colorbar('YTick',(1:no_f_bins)+.5,'YTickLabel',freq_legend)
    
    legend(h1,chan_legend)
    
%     subplot(1,3,3)
%     
%     polar_lim=polar(0, max(max(max(sig_count)), max(abs(sig_vec))), '.'); % Dummy plot to set radial axis length.
% 
%     set(polar_lim,'Marker','none'); % Dummy plot disappears.
%     
%     hold on
%     
%     for ch = 1:2
%         
%         h1(ch) = polar(phase_for_plot, sig_count(:,ch), '-');
%         
%         set(h1(ch),'Color',c_order(ch,:),'LineWidth',2)
%         
%         h = compass(sig_vec(ch));
%         
%         set(h,'Color',c_order(ch,:),'LineWidth',2)
%         
%     end
%     
%     title([folder,', GC Significance Count, by Phase Diff. (S - M)'])
%     
% %     colormap(c_order)
% %     
% %     colorbar('YTick',(1:no_f_bins)+.5,'YTickLabel',freq_legend)
%     
%     legend(h1,chan_legend)
    
    save_as_pdf(gcf,[master_listname(1:end-5),'_GC'])
    
    cd (present_dir)
    
end