function MM_GC_tsplot(filenames, sampling_freq, epoch_length, time_step, freq_lim)

% 'filenames' is a list of filenames containing data previously epoched and 
% run through GC analysis.
% 'sampling_freq' is the sampling frequency of the data.
% 'epoch_length' is the length of each epoch (in seconds, say, for sampling freq. in Hz).
% 'time_step' (optional, can be empty) gives the time step between epochs;
% if time_step is empty, the time step between epochs is epoch_length.
% 'freq_lim' is n by 2, column 1 containing lower cutoff freq., column 2
% containing upper cutoff freq., for n bands over which spectral GC is
% summed.

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

labels = {'Striatal --> Motor','Motor --> Striatal','Difference S->M - M->S'};

color_order = {'k','b','g','m','c','r'};

for file_no = 1:length(filenames)
    
    filename = filenames{file_no};
    
    if isempty(time_step)
        
        listname = [filename,'_channels_',num2str(epoch_length),'s_epochs'];
        
    else
        
        listname = [filename,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epochs'];
        
    end
    
    clear All_GC_spec
    
    GCname = [listname,'_GCspec.mat'];
    
    load(GCname)
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    epoch_numbers = text_read([listname(1:end-1),'_numbers.list'],'%d%*[^\n]');
    
    figure;
    
    %% Plotting by bands.
    
    for a = 1:3
        
        subplot(3,1,a)
    
        for b = 1:no_bands
            
            Bandwise_GC_spec = sum(All_GC_spec(:,f_indices{b},:),2);
            
            t = epoch_numbers*epoch_length;
            
            hold on
            
            h(b) = plot(t,Bandwise_GC_spec(:,:,3),[color_order{b},'-*']);
            
        end
        
        axis tight, box off
        
        title([folder,', Observed Spectral Granger, ',labels{a},', Integrated'])
        
        legend(h,f_lim_label)
        
        plot(t,zeros(size(t)),'k--')
        
        if a == 3
            
            xlabel('Time (s)')
            
        end
        
    end
        
    save_as_pdf(gcf,[GCname(1:end-4),'_',all_f_lim_label,'_ts_diff_only'])
    
end