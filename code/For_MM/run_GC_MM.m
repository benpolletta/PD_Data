function run_GC_MM(filenames, sampling_freq, epoch_length, time_step, maxorder)

% 'filenames' is a list of filenames containing data previously epoched, to
% be run through GC analysis.
% 'sampling_freq' is the sampling frequency of the data.
% 'epoch_length' is the length of each epoch (in seconds, say, for sampling freq. in Hz).
% 'time_step' (optional, can be empty) gives the time step between epochs;
% if time_step is empty, the time step between epochs is epoch_length.
% 'maxorder' is a scalar, determining the maximum order of the VAR. Set it
% high at first - say, 100 for 5000 datapoints; 150 for 20,000; 300 for
% 180,000.

f_length = sampling_freq*epoch_length + 1;

for file_no = 1:length(filenames)
    
    filename = filenames{file_no};
    
    if isempty(time_step)
    
        listname = [filename,'_channels_',num2str(epoch_length),'s_epochs.list'];
    
    else
       
        listname = [filename,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epochs.list'];
        
    end
        
    epoch_list = text_read(listname,'%s%*[^\n]');
    
    no_epochs = length(epoch_list);
    
    All_GC_spec = nan(no_epochs,f_length,2);
    
    Orders = zeros(no_epochs,1);
    Errors = zeros(no_epochs,1);
    
    parfor e = 1:no_epochs
       
        epoch_name = epoch_list{e};
        
        if isempty(dir([epoch_name(1:end-4),'_GC.mat']))
            
            data = load(epoch_name);
            
            [moAIC,info,f] = mvgc_analysis(data',maxorder,epoch_name(1:end-4));
            
        else
            
            mvgc_results = load([epoch_name(1:end-4),'_GC.mat']);
            moAIC = mvgc_results.moAIC;
            info = mvgc_results.info;
            f = mvgc_results.f;
            
        end
        
        Orders(e) = moAIC;
        
        Errors(e) = info.error;
        
        if info.error == 1
            
            display([epoch_name,': ',info.errmsg])
            
        else
            
            GC_spec = reshape([reshape(f(2,1,:),1,f_length) reshape(f(1,2,:),1,f_length)],1,f_length,2);
            
            All_GC_spec(e,:,:) = GC_spec;
            
        end
        
    end
    
    save([listname(1:end-5),'_GCspec.mat'],'Orders','Errors','All_GC_spec')
    
end