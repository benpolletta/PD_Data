function PD_sliding_window_baseline(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin)
    
% Performs sliding window analysis on carbachol data at times of highest
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
%     varargin: further arguments to function_handle.
% SAMPLE CALLS:
% >> load('seven_bands')
% >> data_labels_struct = init_data_labels(freq, no_cycles, bands,
%       'data_field', 'data_subtracted')
% >> subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat',
%       'st_m1_ali2_subjects.mat'};
% Granger causality analysis:
% >> sliding_window_cell = {[500 500], [2 2]};
% >> function_handle = @mvgc_analysis;
% >> varargin = {[], '', 1};
% >> PD_sliding_window_pre_post(function_handle, sliding_window_cell,
%       subjects_mat_cell, data_labels_struct, varargin{:})
% PMTM:
% >> sliding_window_cell = {[500 500], [1 1]};
% >> function_handle = @pmtm;
% >> varargin = {150, [], 500};
% >> PD_sliding_window_pre_post(function_handle, sliding_window_cell,
%       subjects_mat_cell, data_labels_struct, varargin{:})

close('all')

for s = 1:length(subjects_mat_cell)
        
    subjects_mat_name = subjects_mat_cell{s};
    
    load(subjects_mat_name)
    
    no_folders = length(folders);
    
    subjects_mat_prefix = subjects_mat_name(1:(end - length('_subjects.mat')));
    
    parfor fo = 1:no_folders
        
        subjects_struct = load(subjects_mat_name);
        
        folder = subjects_struct.folders{fo};
        
        prefix = subjects_struct.prefixes{fo};
        
        subj_name = [folder, '/', prefix];
        
        data = load([subj_name, get_file_suffix(data_labels_struct, data_labels_struct.data_field)]);
        
        if isfield(data, data_labels_struct.data_field)
            
            data = data.(data_labels_struct.data_field);
            
            display(sprintf('Evaluating %s for %s.', func2str(function_handle), folder))
            
        else
            
            display(sprintf('Data not found for %s.', folder))
            
        end
        
        t = (1:size(data, 1))/data_labels_struct.sampling_freq{1} - subjects_struct.basetimes(fo); % t = t/60;
        
        baseline_data = data(t <= 0, :);
        
        [~, ~] = sliding_window_analysis_multichannel(function_handle, baseline_data,...
            data_labels_struct.sampling_freq, sliding_window_cell, 1, 1, [subj_name, '_baseline',...
            '_band', num2str(data_labels_struct.band_index)], varargin{:});
        
    end
    
end

end