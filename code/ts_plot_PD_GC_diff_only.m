function ts_plot_PD_GC_diff_only(epoch_length, freq_lim)

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
    
    f_lim_label{b} = sprintf('%g-%gHz',freq_lim(b,1),freq_lim(b,2));
    
    all_f_lim_label = [all_f_lim_label,'_',f_lim_label{b}];
    
end

% Grand_GC_name = 'Grand_GC_spec.txt';
%
% Grand_GC_fid = fopen(Grand_GC_name,'w');
%
% GC_format = makeformat(f_length,'f');

color_order = {'k','b','g','m','c','r'};
c_map = [0 0 0; 0 0 1; 0 1 0; 1 1 0; 0 1 1; 1 0 0];

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
    
    cd (folder) 
    
    listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs'];
    
    clear All_GC_spec h
    
    obs_name = [listname,'_GCspec.mat'];
    
    load(obs_name)
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    epoch_numbers = text_read([listname(1:end-1),'_numbers.list'],'%d%*[^\n]');
    
    figure;
    
    for b = 1:no_bands
        
        %% Observed.
        
        Bandwise_GC_spec = sum(All_GC_spec(:,f_indices{b},:),2);
        
        t = epoch_numbers*epoch_length;
        
        hold on
        
        h(b) = plot(t,Bandwise_GC_spec(:,:,3),[color_order{b},'-*']);
        
%         if b==no_bands, colormap(flipud(c_map(1:no_bands,:))), colorbar('YTick',(1:no_bands)+.5,'YTickLabel',fliplr(f_lim_label)), end
        
    end
    
    axis tight, box off
        
    title([folder,', Observed Spectral Granger, S->M - M->S, Integrated'])
    
    legend(h,f_lim_label)
    
    plot(t,zeros(size(t)),'k--')
    
    xlabel('Time (s)')
    
    save_as_pdf(gcf,[obs_name(1:end-4),'_',all_f_lim_label,'_ts_diff_only'])
    
    cd (present_dir)
    
end