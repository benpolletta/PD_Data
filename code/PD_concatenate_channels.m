function PD_concatenate_channels(prefix,plot_opt)

datalist = [prefix,'_datafiles.list'];
    
if isempty(dir(datalist))
    
    file_fid = fopen(datalist,'w');
    
    for digit = 0:9
        
        file = dir([prefix,'*',num2str(digit),'.mat']);
        
        if ~isempty(file)
            
            fprintf(file_fid,'%s\n',file.name);
            
        end
        
    end
    
    fclose(file_fid);
    
end

file_fid = fopen(datalist,'r');

files = textscan(file_fid,'%s');
datafiles = files{1};
datafiles = sort(datafiles);
    
no_files = length(datafiles);
PD_data = [];

sampling_rate = 20*10^3;

last_index = 0;

for f=1:no_files
    
    load(datafiles{f})
    
    if isstruct(data)
        
        data = [data.ch1; data.ch2]';
        
    end
    
    PD_data(last_index+(1:length(data)),:) = data;
    
    clear data
    
    last_index = length(PD_data);
    
end

save([prefix,'_all_channel_data.mat'],'PD_data','sampling_rate')

if plot_opt > 0
   
    figure;
    
    plot((1:length(PD_data))'/(60*sampling_rate),PD_data)
    
    title([prefix,', All Data'])
    
    legend({'Channel 1','Channel 2'})
    
    xlabel('Time (min.)')
    
    save_as_pdf(gcf,[prefix,'_all_channel_data'])
    
end