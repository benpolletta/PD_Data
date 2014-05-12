function plot_PD_GC(epoch_length,freq_lim)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

f_length = sampling_freq*epoch_length + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

if isscalar(freq_lim)

    f_indices = find(f < freq_lim);

    f_lim_label = sprintf('%g-%g',0,freq_lim(1));
    
else
    
    f_indices = find(f > freq_lim(1) & f < freq_lim(2));
    
    f_lim_label = sprintf('%g-%g',freq_lim(1),freq_lim(2));
    
end

% Grand_GC_name = 'Grand_GC_spec.txt';
% 
% Grand_GC_fid = fopen(Grand_GC_name,'w');
% 
% GC_format = makeformat(f_length,'f');

color_order = {'k','b','g','m','c','r'};

labels = {'Striatal --> Motor','Motor --> Striatal','Difference S->M - M->S','Difference z-Scores S->M - M->S'};

for fo = 2:length(folders)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    listname = [prefix,'_channels_',num2str(epoch_length),'s_epochs'];
    
    %% Shuffles.
    
    shuf_name = [listname,'_GCspec_shuffled.mat'];
    
    load(shuf_name)
    
    All_GC_spec = abs(All_GC_spec);
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    % Boxplot, to examine normality.
    
    figure;
    
    for i = 1:3
        
        subplot(3,1,i)
        
        boxplot(All_GC_spec(:,f_indices,i),'plotstyle','compact')
        
        set(gca,'XTick',[1 f_indices(end)],'XTickLabel',f([1 f_indices(end)]))
        
        title([folder,' Spectral Granger, Shuffled, ',labels{i}])
        
    end
    
    xlabel('Frequency (Hz)')
    
    save_as_pdf(gcf,[shuf_name(1:end-4),'_',f_lim_label,'_boxplot'])
    
    % Plotting mean & std.
    
    load([shuf_name(1:end-4),'_stats.mat'])
    
    figure;

    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),shuffled_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),shuffled_mean(:,f_indices,i)+shuffled_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),shuffled_mean(:,f_indices,i)-shuffled_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title([folder,' Spectral Granger, Shuffled, Mean \pm STD'])
    
    save_as_pdf(gcf,[shuf_name(1:end-4),'_',f_lim_label])
    
    %% Observed.
    
    clear All_GC_spec
    
    obs_name = [listname,'_GCspec.mat'];
    
    load(obs_name,'Orders')
    
    All_GC_spec = abs(All_GC_spec);
   
    % Plotting model order.
    
    figure(3)
    
    order_handle(fo) = plot(Orders,['*',color_order{fo}]);
    
    hold on
    
    % Plotting mean, std in freq. domain.
    
    load([obs_name(1:end-4),'_stats.mat'])
    
    figure;

    subplot(2,1,1)
    
    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),spec_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),spec_mean(:,f_indices,i)+spec_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),spec_mean(:,f_indices,i)-spec_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title([folder,' Spectral Granger, Observed, Mean \pm STD'])
        
    subplot(2,1,2)
    
    hold on
    
    for i = 1:4
        
        h(i) = plot(f(f_indices),zs_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),zs_mean(:,f_indices,i)+zs_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),zs_mean(:,f_indices,i)-zs_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title([folder,' Spectral Granger, z-Scored, Mean \pm STD'])
    
    save_as_pdf(gcf,[obs_name(1:end-4),'_',f_lim_label])
    
    cd (present_dir)
    
end

figure(3)

title('VAR Order')
xlabel('Epoch Number')
try
    legend(order_handle,folders)
catch error
    display(error.message)
end
save_as_pdf(gcf,[num2str(epoch_length),'s_epochs_VAR_order'])