function SW_xPlt = PD_pre_post_PLV_xPlt(subjects_mat_cell, filename, norm, data_labels_struct)

if nargin < 4, data_labels_struct = []; end
if isempty(data_labels_struct)
    seven_bands_struct = load('seven_bands');
    data_labels_struct = init_data_labels(seven_bands_struct.freqs, seven_bands_struct.no_cycles, seven_bands_struct.bands, 'data_suffix', '');
end

no_bands = size(data_labels_struct.bands, 1);

short_band_labels = cell(no_bands, 1);

for b = 1:no_bands, short_band_labels{b} = sprintf('%d-%dHz', data_labels_struct.bands(b, :)); end

pd_names = {'pre', 'post'}; 
pd_label = '';
for pd = 1:length(pd_names)
    pd_label = [pd_label, '_', pd_names{pd}];
end

axes_info_struct.names = {'Freq. (Hz)', 'Window_Dim_1', 'Recording', 'Period'};
axes_info_struct.values = {1:200, [], [], pd_names};

for s = 1:length(subjects_mat_cell)
    
    subjects_mat_struct = load(subjects_mat_cell{s});
    
    axes_info_struct.values{strcmp(axes_info_struct.names, 'Recording')} = subjects_mat_struct.folders;

    subjects_mat_name = subjects_mat_cell{s}(1:(end - length('_subjects.mat')));
    
    load([subjects_mat_name, data_labels_struct.BP_suffix, data_labels_struct.peak_suffix,...
        '_pct_', short_band_labels{data_labels_struct.band_index}, '_high_', num2str(data_labels_struct.epoch_secs/60),...
        '_min_secs', data_labels_struct.pd_handle, '_PLV.mat'])

    sw_size = size(Coh_sec_pct);
    
    sw_xPlt = xPlt;

    sw_xPlt = sw_xPlt.importData({Coh_sec_pct});

    for dim = length(axes_info_struct.names):-1:2
        
        if sw_size(dim) < 100
        
            sw_xPlt = unpackDim(sw_xPlt, dim, 1, axes_info_struct.names{dim}, axes_info_struct.values{dim});
        
        end
            
    end
    
    if exist('SW_xPlt', 'var')
        
        SW_xPlt = SW_xPlt.merge(sw_xPlt);
        
    else
        
        SW_xPlt = sw_xPlt;
        
    end
            
end

for dim = length(axes_info_struct.names):-1:2
    
    if sw_size(dim) >= 100
        
        SW_xPlt = unpackDim(SW_xPlt, dim, dim, axes_info_struct.names{dim}, axes_info_struct.values{dim});
        
    end
    
end

SW_xPlt.axis((end - 1):end) = [];

SW_xPlt.meta.matrix_dim_1 = nDDictAxis;

SW_xPlt.meta.matrix_dim_1.name = axes_info_struct.names{1};
SW_xPlt.meta.matrix_dim_1.values = axes_info_struct.values{1};

save([filename, pd_label,'_band', num2str(data_labels_struct.band_index),...
    data_labels_struct.BP_suffix, data_labels_struct.peak_suffix, data_labels_struct.pd_handle, norm, '_PLV'],...
    'SW_xPlt', 'axes_info_struct')
    
end