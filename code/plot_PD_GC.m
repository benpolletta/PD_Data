function plot_PD_GC

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

epoch_length = 5*sampling_freq;

f_length = epoch_length/2 + 1;

f = (sampling_freq/2)*(1:f_length)/f_length;

f_indices = find(f < 200);

% Grand_GC_name = 'Grand_GC_spec.txt';
% 
% Grand_GC_fid = fopen(Grand_GC_name,'w');
% 
% GC_format = makeformat(f_length,'f');

color_order = {'k','b','g','m','c','r'};

labels = {'Striatal --> Motor','Motor --> Striatal','Difference S->M - M->S'};

for fo = 2:(length(folders)-1)
    
    folder = folders{fo};
    
    prefix = prefixes{fo};
   
    cd (folder)
    
    %% Shuffles.
    
    shuf_name = [prefix,'_channels_epochs_GCspec_shuffled.mat'];
    
    load(shuf_name)
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    % Boxplot, to examine normality.
    
    figure;
    
    for i = 1:3
        
        subplot(3,1,i)
        
        boxplot(All_GC_spec(:,f_indices,i),'plotstyle','compact')
        
        set(gca,'XTick',[1 f_indices(end)],'XTickLabel',f([1 f_indices(end)]))
        
        title(['Spectral Granger, Shuffled, ',labels{i}])
        
    end
    
    xlabel('Frequency (Hz)')
    
    save_as_pdf(gcf,[shuf_name(1:end-4),'_boxplot'])
    
    % Calculating, plotting mean & std.
    
    shuffled_mean = mean(All_GC_spec);
    
    shuffled_std = std(All_GC_spec);
    
    figure;

    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),shuffled_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),shuffled_mean(:,f_indices,i)+shuffled_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),shuffled_mean(:,f_indices,i)-shuffled_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title('Spectral Granger, Shuffled, Mean \pm STD')
    
    save_as_pdf(gcf,shuf_name(1:end-4))
    
    %% Observed.
    
    obs_name = [prefix,'_channels_epochs_all_GCspec.mat'];
    
    load(obs_name)
   
    % Plotting model order.
    
    figure(3)
    
    order_handle(fo) = plot(Orders,['*',color_order{fo}]);
    
    hold on
    
    % Computing difference of GC.
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
    spec_mean = mean(All_GC_spec);
    
    spec_std = std(All_GC_spec);
    
    % Computing normalized GC.
    
    All_GC_norm = nan(size(All_GC_spec));
    
    for i=1:3
        
        All_GC_norm(:,:,i) = (All_GC_spec(:,:,i) - ones(size(All_GC_spec(:,:,i)))*diag(shuffled_mean(:,:,i)))*diag(shuffled_std(:,:,i));
        
    end
    
    save([obs_name(1:end-4),'_thresh.mat'],'All_GC_norm')
    
    norm_mean = mean(All_GC_norm);
    
    norm_std = std(All_GC_norm);
    
    figure;

    subplot(2,1,1)
    
    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),spec_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),spec_mean(:,f_indices,i)+spec_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),spec_mean(:,f_indices,i)-spec_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title('Spectral Granger, Observed, Mean \pm STD')
        
    subplot(2,1,2)
    
    hold on
    
    for i = 1:3
        
        h(i) = plot(f(f_indices),norm_mean(:,f_indices,i),color_order{i});
        
        plot(f(f_indices),norm_mean(:,f_indices,i)+norm_std(:,f_indices,i),[':',color_order{i}])
        
        plot(f(f_indices),norm_mean(:,f_indices,i)-norm_std(:,f_indices,i),[':',color_order{i}])
        
    end
    
    legend(h,labels)
    
    title('Spectral Granger, z-Scored, Mean \pm STD')
    
    save_as_pdf(gcf,obs_name(1:end-4))
    
    cd (present_dir)
    
end

figure(3)

title('VAR Order')
xlabel('Epoch Number')
legend(order_handle,folders)