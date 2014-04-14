present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

for fo = 1:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    for pd = -1:2
    
        filename = dir([folder,'/',prefix,'_5min',num2str(pd),'/AVG_THRESH_FILE_SHUFFLE_',prefix,'_5min',num2str(pd),'_1000_shufs/*IE_zs.fig']);
        
        if ~isempty(filename)
        
            open([folder,'/',prefix,'_5min',num2str(pd),'/AVG_THRESH_FILE_SHUFFLE_',prefix,'_5min',num2str(pd),'_1000_shufs/',filename.name])
        
        end
            
    end
    
end

titles = {'Pre1','Post1','Post2','Post3'};
xlabels = {'Phase (Hz)','Phase (Hz)','Phase (Hz)','Phase (Hz)'};
ylabels = {'Mouse 1','Mouse 2','Mouse 3','Mouse 4','Mouse 5'};

[max_data_all,data_all] = figure_replotter_labels(1:20,5,4,5,5,1:.25:30,20:5:250,titles,xlabels,ylabels);

save_as_pdf(gcf,'PD_all_5min_thresh')

figure;

for i=1:4
    
    median_data_all = nanmedian(data_all(:,:,(0:4:19)+i),3);
    
    subplot(1,4,i)
    
    imagesc(1:.25:30,20:5:250,median_data_all)
    
    axis xy
    
    xlabel('Phase (Hz)')
    ylabel('Amplitude (Hz)')
    title([num2str((i-2)*5),' to ',num2str((i-2)*5 + 4),' Min. Rel. Inj.'])
    
end

save_as_pdf(gcf,'PD_5min_thresh_median')