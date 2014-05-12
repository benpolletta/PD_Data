function MM_epochs(filenames, sampling_freq, epoch_length, time_step, outlier_check)

% 'filenames' is a list of filenames containing data to be epoched.
% 'sampling_freq' is the sampling frequency of the data.
% 'epoch_length' is the length of each epoch (in seconds, say, for sampling freq. in Hz).
% 'time_step' (optional, can be empty) gives the time step between epochs;
% if time_step is empty, the time step between epochs is epoch_length.
% 'outlier_check' is the max. standard deviation (relative to the entire
% data segment) that is allowed before an epoch is considered an outlier.

for file_no = 1:length(filenames)
    
    filename = filenames{file_no};
    
    all_data = load(filename);
    
    [r,c] = size(all_data);
    
    if r < c
        
        all_data = all_data';
        
    end
    
    [data_length, no_channels] = size(all_data);
    
    data_format = makeformat(no_channels,'f');
    
    if isempty(time_step)
        
        listname = [filename,'_channels_',num2str(epoch_length),'s'];
    
        no_epochs = floor(data_length/(epoch_length*sampling_freq));
        
        time_step = epoch_length;
        
    else
        
        listname = [filename,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s'];
    
        no_epochs = floor((data_length - epoch_length*sampling_freq)/(time_step*sampling_freq));
        
    end
    
    %% Set up for artifact flagging.
    
    all_mean = mean(all_data);
    
    all_std = std(all_data);
    
    %% Open files for lists.
    
    outliers_fid = fopen([listname,'_outliers.list'],'w');
    
    epoch_list_fid = fopen([listname,'_epochs.list'],'w');
    
    epoch_nos_list_fid = fopen([listname,'_epoch_numbers.list'],'w');
    
    %% Extracting Epochs.
    
    for e = 1:no_epochs
        
        epoch_name = [listname,'_epoch',num2str(e),'.txt'];
        
        % Getting epoch.
        
        start_sec = (e-1)*time_step;
        
        end_sec = (e-1)*time_step + epoch_length;
        
        start_index = max(start_sec*sampling_freq + 1, 1);
        
        end_index = min(end_sec*sampling_freq, length(all_data));
        
        if start_index < end_index
            
            data = all_data(start_index:end_index,:);
            
            fid = fopen(epoch_name, 'w');
            
            fprintf(fid, data_format, data');
            
            fclose(fid);
            
            % Determining whether data contains outliers.
            
            norm_data = (data - ones(size(data))*diag(all_mean))*diag(1./all_std);
            
            if any(abs(norm_data) > outlier_check)
                
                zs_name = [epoch_name(1:end-4),'_',num2str(epoch_length),'s_by_',num2str(time_step),'s_zs.txt'];
                
                fid = fopen(zs_name, 'w');
                
                fprintf(fid, data_format, norm_data');
                
                fclose(fid);
                
                fprintf(outliers_fid,'%s\n',zs_name);
                
            else
                
                fprintf(epoch_nos_list_fid,'%d\n',e);
                
                fprintf(epoch_list_fid,'%s\n',epoch_name);
                
            end
            
        end
        
    end

    fclose('all');
    
end