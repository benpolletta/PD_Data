function MM_GC_colorplot(filenames, sampling_freq, epoch_length, time_step, freq_lim)

epoch_datapoints = epoch_length*sampling_freq;

f_length = epoch_datapoints + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

if isscalar(freq_lim)

    f_indices = find(f < freq_lim);

    f_lim_label = sprintf('%g-%g',0,freq_lim(1));
    
else
    
    f_indices = find(f > freq_lim(1) & f < freq_lim(2));
    
    f_lim_label = sprintf('%g-%g',freq_lim(1),freq_lim(2));
    
end

labels = {'Striatal --> Motor','Motor --> Striatal','Difference S->M - M->S'};

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
    
    %% Plotting parameters regarding VAR fit.
    
    figure;
    
    subplot(2,1,1)
    
    plot(Orders,'*')
    
    title(listname)
    
    ylabel('VAR Order')
    
    subplot(2,1,2)
    
    plot(Orders,'*')
    
    xlabel('Epoch Number')
    
    ylabel('Error in VAR Fitting?')
    
    save_as_pdf(gcf,[GCname(1:end-4),'_fit'])
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    %% Plotting spectral GC.
    
    figure;
    
    for i = 1:3
        
        subplot(3,1,i)
        
        imagesc((1:size(All_GC_spec,1))*epoch_length,f(f_indices),All_GC_spec(:,f_indices,i)')
        
%         caxis([median(min(All_GC_spec(:,f_indices,i),[],2)) median(max(All_GC_spec(:,f_indices,i),[],2))])
        
        ylabel('Frequency (Hz)')
        
        axis xy, colorbar
        
        title([folder,' Spectral Granger, Observed, ',labels{i}])
        
    end
    
    xlabel('Time (s)')
    
    save_as_pdf(gcf,[GCname(1:end-4),'_',f_lim_label,'_cplot'])
    
end