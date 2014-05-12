function colorplot_PD_GC_noshuffle(epoch_length,freq_lim)

present_dir = pwd;

folders = {'090312','130703','130709','130718','130725'};

prefixes = {'09312','13703','13709','13718','13725'};

sampling_freq = 1000;

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
    
% Grand_GC_name = 'Grand_GC_spec.txt';
% 
% Grand_GC_fid = fopen(Grand_GC_name,'w');
% 
% GC_format = makeformat(f_length,'f');

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
    
    All_GC_spec(:,:,3) = -diff(All_GC_spec,[],3);
    
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
    
    save_as_pdf(gcf,[obs_name(1:end-4),'_',f_lim_label,'_cplot'])
    
    cd (present_dir)
    
end