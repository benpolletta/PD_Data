function PD_concatenate_channels(prefix, plot_opt, file_format)

datalist = [prefix,'_datafiles.list'];
    
% if isempty(dir(datalist))
    
    file_fid = fopen(datalist,'w');
    
    for digit = 0:9
    
        if strcmp(file_format, 'mat') || isempty(file_format)
            
            file = dir([prefix, '*', num2str(digit), '.mat']);
            
        elseif strcmp(file_format(1:3), 'abf')
            
            file = dir([prefix, '*', num2str(digit), '.abf']);
        
        end
            
        if ~isempty(file)
            
            fprintf(file_fid, '%s\n', file.name);
            
        end
        
    end
    
    fclose(file_fid);
    
% end

file_fid = fopen(datalist,'r');

files = textscan(file_fid,'%s');
datafiles = files{1};
datafiles = sort(datafiles);
    
no_files = length(datafiles);
PD_data = [];

sampling_rate = 20*10^3;

last_index = 0;

for f=1:no_files
    
    if strcmp(file_format, 'mat') || isempty(file_format)
        
        load(datafiles{f})
        
        if isstruct(data)
            
            data = [data.ch1; data.ch2]';
            
        end
    
    elseif strcmp(file_format, 'abf')
        
        data = abf2load(datafiles{f});
        
    elseif strcmp(file_format, 'abf_ali')
        
        data = abf2load(datafiles{f});
        
        data = permute(data(:, 1:2, :), [1 3 2]);
        
        data_size = size(data);
        
        data = reshape(data, data_size(1)*data_size(2), data_size(3));
        
    end
    
    PD_data(last_index + (1:length(data)), :) = data;
    
    clear data
    
    last_index = length(PD_data);
    
end

save([prefix,'_all_channel_data.mat'],'-v7.3','PD_data','sampling_rate')

if plot_opt > 0
   
    figure;
    
    plot((1:length(PD_data))'/(60*sampling_rate), PD_data)
    
    title([prefix,', All Data'])
    
    legend({'Channel 1','Channel 2'})
    
    xlabel('Time (min.)')
    
    save_as_pdf(gcf,'-v7.3',[prefix,'_all_channel_data'])
    
end