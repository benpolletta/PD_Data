%% General stuff

load('seven_bands')

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted')

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

filename = 'STR_w_M1';

%% Granger causality analysis:

sliding_window_cell = {10*500*[1 1], [2 2]};

% window_size = cellfun(@(x,y) x./y, sliding_window_cell, data_labels_struct.sampling_freq);

function_handle = @mvgc_analysis; function_name = function_handle;

varargin = {[], '', 1};

% %% PMTM:
% 
% sliding_window_cell = {[150*500 150*500], [1 1]};
% 
% function_handle = @pmtm;
% 
% varargin = {150, [], 500};

% %% AR spectrum:
% 
% sliding_window_cell{2} = [1 1]; % {[500 500], [1 1]}; % 
% 
% function_handle = @mvgc_analysis;
% 
% varargin = {[], '', 3};

% %% AR cross-spectrum:
% 
% sliding_window_cell = {[150*500 150*500], [2 2]};
% 
% function_handle = @mvgc_analysis;
% 
% varargin = {[], '', 3};

PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

% total_dims = sum(window_size > 1) + sum(output_size > 1);
% 
% channel_dims = [2 3];
% 
% indices_cell = make_indices_cell(total_dims, channel_dims, subjects_mat_cell);

% SW = PD_sliding_window_pre_post_load(function_name, sliding_window_cell, subjects_mat_cell, data_labels_struct, indices_cell, filename, varargin{:});