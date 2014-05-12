function ts_plot_PD_GC(epoch_length, freq_lim)

% freq_lim is n by 2, column 1 containing lower cutoff freq., column 2
% containing upper cutoff freq., for n bands over which spectral GC is
% summed.

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

epoch_datapoints = epoch_length*sampling_freq;

f_length = epoch_datapoints + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

no_bands = size(freq_lim,1);

all_f_lim_label = '';

for b = 1:no_bands
    
    f_indices{b} = find(f > freq_lim(b,1) & f < freq_lim(b,2));
    
    f_lim_label{b} = sprintf('%g-%g',freq_lim(b,1),freq_lim(b,2));
    
    all_f_lim_label = [all_f_lim_label,'_',f_lim_label{b}];
    
end

% Grand_GC_name = 'Grand_GC_spec.txt';
%
% Grand_GC_fid = fopen(Grand_GC_name,'w');
%
% GC_format = makeformat(f_length,'f');

labels = {'Striatal --> Motor','Motor --> Striatal','S->M - M->S','S->M z-sc. - M->S z-sc.'};

color_order = {'k','b','g','m','c','r'};
c_map = [0 0 0; 0 0 1; 0 1 0; 1 1 0; 0 1 1; 1 0 0];

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    cd (folder)    
    
    listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs'];
    
    epoch_numbers = text_read([listname(1:end-1),'_numbers.list'],'%d%*[^\n]');
    
    figure;
    
    for b = 1:no_bands
        
        %% Shuffles.
        
        shuf_name = [listname,'_GCspec_shuffled.mat'];
        
        load(shuf_name)
        
        All_GC_spec = abs(All_GC_spec);
        
        All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
        
        Bandwise_GC_spec = sum(All_GC_spec(:,f_indices{b},:),2);
        
        subplot(3,no_bands,b)
        
        hold on
        
        clear h
        
        for i = 1:3
            
            h(i) = plot(Bandwise_GC_spec(:,:,i),[color_order{i},'*']);
            
        end
        
        % legend(h,labels,'location','NorthEastOutside')
        
        if b==no_bands, colormap(flipud(c_map(1:3,:))), colorbar('YTick',(1:3)+.5,'YTickLabel',fliplr(labels(1:3))), end
        
        axis tight, box off
        
        title([folder,' Shuffled Spectral Granger, Integrated Over ',f_lim_label{b},' Hz'])
        
        % save_as_pdf(gcf,[shuf_name(1:end-4),'_',f_lim_label{b},'_ts'])
        
        %% Observed.
        
        clear All_GC_spec h
    
        obs_name = [listname,'_GCspec.mat'];
        
        load(obs_name)
        
        All_GC_spec = abs(All_GC_spec);
        
        All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
        
        Bandwise_GC_spec = sum(All_GC_spec(:,f_indices{b},:),2);
        
        t = epoch_numbers*epoch_length;
        
        subplot(3,no_bands,no_bands+b)
        
        hold on
        
        clear h
        
        for i = 1:3
            
            h(i) = plot(t,Bandwise_GC_spec(:,:,i),[color_order{i},'-*']);
            
        end
        
        % legend(h,labels,'location','NorthEastOutside')
        
        if b==no_bands, colormap(flipud(c_map(1:3,:))), colorbar('YTick',(1:3)+.5,'YTickLabel',fliplr(labels(1:3))), end
        
        axis tight, box off
        
        title([folder,', Observed Spectral Granger, Integrated Over ',f_lim_label{b},' Hz'])
        
        % save_as_pdf(gcf,[obs_name(1:end-4),'_',f_lim_label{b},'_ts'])
        
        load([obs_name(1:end-4),'_thresh.mat'])
        
        Bandwise_GC_zs = sum(All_GC_zs(:,f_indices{b},:),2);
        
        t = epoch_numbers*epoch_length;
        
        subplot(3,no_bands,2*no_bands+b)
        
        hold on
        
        for i = 1:3
            
            h(i) = plot(t,Bandwise_GC_zs(:,:,i),[color_order{i},'-*']);
            
        end
        
        % legend(h,labels,'location','NorthEastOutside')
        
        if b==no_bands, colormap(flipud(c_map(1:3,:))), colorbar('YTick',(1:3)+.5,'YTickLabel',fliplr(labels(1:3))), end
        
        axis tight, box off
        
        title([folder,', z-Scored Spectral Granger, Integrated Over ',f_lim_label{b},' Hz'])
        
        save([obs_name(1:end-4),'_',f_lim_label{b},'_ts'],'Bandwise_GC_spec','Bandwise_GC_zs')
        
    end
    
    save_as_pdf(gcf,[obs_name(1:end-4),'_',all_f_lim_label,'_ts'])
    
    cd (present_dir)
    
end