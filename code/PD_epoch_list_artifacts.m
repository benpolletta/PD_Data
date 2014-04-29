function PD_epoch_list_artifacts(prefix, injection_sec, outlier_check)

load([prefix,'_all_data_dec.mat'])

%% Set up for epoching.

t_secs = 1:length(PD_dec)/sampling_freq;

epoch_length = 6;

no_pre = ceil(injection_sec/epoch_length);
no_post = ceil(t_secs(end)/epoch_length);

%% Set up for artifact flagging.

all_mean = mean(PD_dec);
all_std = std(PD_dec);

%% Set up for lists.

outliers_fid = fopen([prefix,'_outliers.list'],'w');

master_fid = fopen([prefix,'_min_master.list'],'w');
min_epochs_fid = fopen([prefix,'_min_epochs.list'],'w');

min_list_fid = nan(length(-5:14),1);
min_lims = nan(length(-5:14),2);

for l = -5:14
    
    listname = [prefix,'_min',num2str(l),'.list'];
    
    min_list_fid(l+6) = fopen(listname,'w');
    fprintf(master_fid,'%s\n',listname);
    
    min_lims(l+6,:) = [l l+1];
    
end

fclose(master_fid);

master_fid = fopen([prefix,'_5min_master.list'],'w');
five_min_epochs_fid = fopen([prefix,'_5min_epochs.list'],'w');

five_min_list_fid = nan(length(-1:2),1);
five_min_lims = nan(length(-1:2),2);

for l = -1:2
    
    listname = [prefix,'_5min',num2str(l),'.list'];
    
    five_min_list_fid(l+2) = fopen(listname,'w');
    fprintf(master_fid,'%s\n',listname);
    
    five_min_lims(l+2,:) = [l*5 (l+1)*5];
    
end

fclose(master_fid);

%% Extracting Epochs.

for e = -no_pre:(no_post-1)
    
    epoch_name = [prefix,'_epoch',num2str(e),'.txt'];
    
    % Getting epoch.
   
    start_sec = injection_sec + e*epoch_length;
    
    end_sec = injection_sec + (e+1)*epoch_length;
    
    start_index = max(start_sec*sampling_freq + 1, 1);
    
    end_index = min(end_sec*sampling_freq, length(PD_dec));
    
    if start_index < end_index
        
        data = PD_dec(start_index:end_index);
        
        fid = fopen(epoch_name,'w');
        
        fprintf(fid,'%f\n',data);
        
        fclose(fid);
        
        % Determining list indices.
        
        start_min = e*epoch_length/60;
        
        end_min = e*epoch_length/60;
        
        min_index = min_lims(:,1) <= start_min & end_min < min_lims(:,2);
        
        five_min_index = five_min_lims(:,1) <= start_min & end_min < five_min_lims(:,2);
        
        % Determining whether data contains outliers.
        
        norm_data = (data - all_mean)/all_std;
        
        if any(abs(norm_data) > outlier_check)
            
            zs_name = [epoch_name(1:end-4),'_zs.txt'];
            
            fid = fopen(zs_name,'w');
            
            fprintf(fid,'%f\n',norm_data);
            
            fclose(fid)
            
            fprintf(outliers_fid,'%s\n',zs_name);
            
        else
            
            if any(min_index > 0)
                
                fprintf(min_list_fid(min_index),'%s\n',epoch_name);
                
            end
            
            if any(five_min_index > 0)
                
                fprintf(five_min_list_fid(five_min_index),'%s\n',epoch_name);
                
            end
            
        end
        
    end
    
end

fclose('all');