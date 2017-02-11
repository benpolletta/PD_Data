function file_suffix = get_file_suffix(data_labels_struct, data_type)

file_suffix = '';

switch data_type
    
    case 'PD_data'
        
        file_suffix = '_all_channel_data.mat';
        
    case 'PD_dec'
        
        file_suffix = ['_all_channel_data_dec', data_labels_struct.data_suffix, '.mat'];
        
    case 'high_density_periods'
        
        file_suffix = [data_labels_struct.BP_suffix, data_labels_struct.peak_suffix,...
        '_pct_BP_high_', num2str(data_labels_struct.epoch_secs/60), '_min_secs', data_labels_struct.pd_handle, '.mat'];
        
end

if isempty(file_suffix)
    
    display(sprintf('Data type %s has no defined file suffix in get_file_suffix.m.', data_type))
    
end