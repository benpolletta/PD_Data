function plot_PD_GC_noshuffle(epoch_length,freq_lim)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

f_length = sampling_freq*epoch_length/2 + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

if isscalar(freq_lim)

    f_indices = find(f < freq_lim);

    f_lim_label = sprintf('%g-%g',0,freq_lim(1));
    
else
    
    f_indices = find(f > freq_lim(1) & f < freq_lim(2));
    
    f_lim_label = sprintf('%g-%g',freq_lim(1),freq_lim(2));
    
end

color_order = {'k','b','g','m','c','r'};

labels = {'Striatal --> Motor','Motor --> Striatal','Difference S->M - M->S','Difference z-Scores S->M - M->S'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs'];
    
    %% Observed.
    
    clear All_GC_spec
    
    obs_name = [listname,'_GCspec.mat'];
    
    load(obs_name)
    
    All_GC_spec = abs(All_GC_spec);
   
    % Plotting model order.
    
    figure(1)
    
    order_handle(fo) = plot(Orders,['*',color_order{fo}]);
    
    hold on
    
    % Computing difference of spectral GC, mean & median.
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    spec_mean = nanmean(All_GC_spec);
    
    spec_std = nanstd(All_GC_spec);
    
    % Plotting spectral GC.
    
    figure;
    
    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),spec_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),spec_mean(:,f_indices,i)+spec_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),spec_mean(:,f_indices,i)-spec_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title([folder,' Spectral Granger, Observed, Mean \pm STD'])
    
    save_as_pdf(gcf,[obs_name(1:end-4),'_',f_lim_label,'_nothresh'])
    
    cd (present_dir)
    
end

figure(1)

title('VAR Order')
xlabel('Epoch Number')
try
    legend(order_handle,folders)
catch error
    display(error.message)
end
save_as_pdf(gcf,[num2str(epoch_length),'s_epochs_VAR_order'])