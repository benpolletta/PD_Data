function SW_xPlt = PD_sliding_window_load(function_name, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, pd_names, varargin)
    
% Loads sliding window analysis on carbachol data at times of highest
% striatal beta band density, pre- and post-infusion.
% INPUTS:
%     function_handle (function handle): analysis to be performed on each window.
%     sliding_window_cell (2 x 1 cell of 2 x 1 arrays): cell
%       containing sliding window length and step length (in indices) for the
%       two dimensions of the carbachol data.
%     subjects_mat_cell (n x 1 cell of strings): cell containing names of
%       *subjects.mat files to be (batch) analyzed.
%     data_labels_struct (structure): can be initialized by
%       init_data_labels.m, contains fields: BP_suffix, peak_suffix,
%       data_suffix, epoch_secs, pd_handle, data_field, band_index,
%       sampling_freq.
%     indices_cell (n x 2 cell of cells): allows permutation of dimension,
%       for each set of analyses (i.e., each *subjects.mat file) - to
%       handle, e.g., different location of striatal channel for different
%       sets of recordings:
%       >> indices_cell = {{':', [1 2]}, {':', [1 2]};...
%               {':', [1 2]}, {':', [2 1]};...
%               {':', [1 2]}, {':', [2 1]}};       
%     filename (string): name of collected analysis (e.g., 'STR_w_M1').

function_name = get_fname(function_name);
    
window_time_cell = cellfun(@(x,y) x./y, sliding_window_cell, data_labels_struct.sampling_freq, 'UniformOutput', 0);

% pd_names = {'pre', 'post'};

no_periods = length(pd_names);

pd_label = '';

for period = 1:no_periods
    
    pd_label = [pd_label, '_', pd_names{period}];
    
end

for s = 1:length(subjects_mat_cell)
    
    load(subjects_mat_cell{s})
    
    subject_mat_folders(s) = length(folders);
    
end

no_folders = sum(subject_mat_folders);
            
first_name = make_sliding_window_analysis_name([folders{1}, '/', prefixes{1}, '_', pd_names{1},...
    '_band', num2str(data_labels_struct.band_index)], function_name, window_time_cell, 2, varargin{:});

load(first_name, 'output_size', 'no_windows')

indices_cell = make_PD_indices_cell(subjects_mat_cell,function_name,...
    sliding_window_cell, no_windows, output_size, varargin{:});

sw_size = [output_size, no_windows'];

sw_size(sw_size == 1) = [];

sw_size_cell = num2cell(sw_size);

SW = nan(sw_size_cell{:}, no_folders, no_periods);

folder_index = 0;

for s = 1:length(subjects_mat_cell)
    
    subjects_struct = load(subjects_mat_cell{s});
    
    for fo = 1:length(subjects_struct.folders)
        
        folder_index = folder_index + 1;
        
        folder = subjects_struct.folders{fo};
        
        prefix = subjects_struct.prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        for period = 1:no_periods
            
            sw_xPlt = PD_sliding_window_xPlt(function_name, sliding_window_cell,...
                subjects_struct, data_labels_struct, subj_name, pd_names(period), varargin{:});
            
            sw_xPlt.axis(end + 1).name = 'Recording';
            sw_xPlt.axis(end).values = {folder};
            
            sw_xPlt.axis(end + 1).name = 'Period';
            sw_xPlt.axis(end).values = {pd_names{period}};
            
            if exist('SW_xPlt', 'var')
            
                SW_xPlt = SW_xPlt.merge(sw_xPlt);
                
            else
                
                SW_xPlt = sw_xPlt;
                
            end 
            
        end
        
    end
    
end

SW_xPlt = squeeze(SW_xPlt);

save([make_sliding_window_analysis_name([filename, pd_label,...
    '_band', num2str(data_labels_struct.band_index)], function_name,...
    window_time_cell, 2, varargin{:}), '_xPlt.mat'], 'SW_xPlt')