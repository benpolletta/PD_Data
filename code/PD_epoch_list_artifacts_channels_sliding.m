function PD_epoch_list_artifacts_channels_sliding(prefix, epoch_length, time_step, outlier_check)

load([prefix,'_all_channel_data_dec.mat'])

%% Set up for epoching.

data_length = length(PD_dec);

% no_epochs = floor(data_length/(epoch_length*sampling_freq));
no_epochs = floor((data_length - epoch_length*sampling_freq)/(time_step*sampling_freq));

%% Set up for artifact flagging.

all_mean = mean(PD_dec);
all_std = std(PD_dec);

%% Set up for lists.

outliers_fid = fopen([prefix,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_outliers.list'],'w');

epoch_list_fid = fopen([prefix,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epochs.list'],'w');

epoch_nos_list_fid = fopen([prefix,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epoch_numbers.list'],'w');

%% Extracting Epochs.

for e = 1:no_epochs
    
    epoch_name = [prefix,'_channels_',num2str(epoch_length),'s_by_',num2str(time_step),'s_epoch',num2str(e),'.txt'];
    
    % Getting epoch.
   
    start_sec = (e-1)*time_step;
    
    end_sec = (e-1)*time_step + epoch_length;
    
    start_index = max(start_sec*sampling_freq + 1, 1);
    
    end_index = min(end_sec*sampling_freq, length(PD_dec));
    
    if start_index < end_index
        
        data = PD_dec(start_index:end_index,:);
        
        fid = fopen(epoch_name,'w');
        
        fprintf(fid,'%f\t%f\n',data');
        
        fclose(fid);
        
        % Determining whether data contains outliers.
        
        norm_data = (data - ones(size(data))*diag(all_mean))*diag(1./all_std);
        
        if any(abs(norm_data) > outlier_check)
            
            zs_name = [epoch_name(1:end-4),'_',num2str(epoch_length),'s_by_',num2str(time_step),'s_zs.txt'];
            
            fid = fopen(zs_name,'w');
            
            fprintf(fid,'%f\t%f\n',norm_data');
            
            fclose(fid);
            
            fprintf(outliers_fid,'%s\n',zs_name);
            
        else
        
            fprintf(epoch_nos_list_fid,'%d\n',e);
            
            fprintf(epoch_list_fid,'%s\n',epoch_name);
            
        end
        
    end
    
end

fclose('all');