function PD_sliding_window_shuffle(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, pd_names, shuffle_struct, varargin)
    
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

if isempty(shuffle_struct)
    
    shuffle_struct = init_struct({'no_shuffles', 'combine_pairs', 'varargin'},...
        {1000, @combine_dimension, {2, [1 2]}});
    
end

no_pds = length(pd_names);

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
        
        max_beta_density = load([subjects_mat_prefix, get_file_suffix(data_labels_struct, 'high_density_periods')]);
        
        striatal_channel = find(strcmp(subjects_struct.chan_labels, 'Striatum'));
        
        if isfield(data, data_labels_struct.data_field)
            
            data = data.(data_labels_struct.data_field);
            
            display(sprintf('Evaluating %s for %s.', func2str(function_handle), folder))
            
        else
            
            display(sprintf('Data not found for %s.', folder))
            
        end
        
        t = (1:size(data, 1))/data_labels_struct.sampling_freq{1} - subjects_struct.basetimes(fo); % t = t/60;
        
        for pd = 1:no_pds
            
            switch pd_names{pd}
            
                case 'pre'
                    
                    bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, data_labels_struct.band_index, 1);
                    
                    bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, data_labels_struct.band_index, 1);
                    
                    data_index = t >= t(bp_max_start) & t <= t(bp_max_end);
                    
                case 'post'
                    
                    bp_max_start = max_beta_density.All_bp_max_start(fo, striatal_channel, data_labels_struct.band_index, 2);
                    
                    bp_max_end = max_beta_density.All_bp_max_end(fo, striatal_channel, data_labels_struct.band_index, 2);
                    
                    data_index = t >= t(bp_max_start) & t <= t(bp_max_end);

                case 'baseline'

                    data_index = t <= 0;

            end
            
            data_selected = data(data_index, :);
            
            [~, ~] = sliding_window_shuffle_multichannel(function_handle, data_selected,...
                data_labels_struct.sampling_freq, sliding_window_cell, shuffle_struct, 1,...
                [subj_name, '_', pd_names{pd}, '_shuffles_band', num2str(data_labels_struct.band_index)],...
                varargin{:});
            
        end
        
    end
    
end

end