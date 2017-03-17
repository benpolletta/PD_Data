%% General stuff

load('seven_bands')

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted')

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

chan_labels = {'Striatum', 'M1'};

filename = 'STR_w_M1';

window_length = 10;

sliding_window_cell{1} = window_length*500*[1 1];

% %% Granger causality analysis:
% 
% sliding_window_cell{2} = [2 2];
% 
% function_handle = @mvgc_analysis; function_name = function_handle;
% 
% varargin = {[], '', 1};
% 
% % output_struct.unpack_flag = 1;
% % 
% % output_struct.output_names = {'Freq. (Hz)', 'To', 'From'};
% % 
% % output_struct.output_values = {250*(1:(window_length*500 + 1))/(window_length*500), chan_labels, chan_labels};

% %% PMTM:
% 
% sliding_window_cell{2} = [1 1];
% 
% function_handle = @pmtm;
% 
% varargin = {150, [], 500};

%% AR spectrum:

sliding_window_cell{2} = [1 1]; % {[500 500], [1 1]}; % 

function_handle = @mvgc_analysis;

varargin = {[], '', 3};

% output_struct.unpack_flag = 1;
% 
% output_struct.output_names = {'Freq. (Hz)', 'Channel'};
% 
% output_struct.output_values = {250*(1:(window_length*500 + 1))/(window_length*500), chan_labels};

% %% AR cross-spectrum:
% 
% sliding_window_cell{2} = [2 2]; % = {[150*500 150*500], [2 2]};
% 
% function_handle = @mvgc_analysis;
% 
% varargin = {[], '', 3};
% 
% % output_struct.unpack_flag = 1;
% % 
% % output_struct.output_names = {'Freq. (Hz)', 'Channel 1', 'Channel 2'};
% % 
% % output_struct.output_values = {250*(1:(window_length*500 + 1))/(window_length*500), chan_labels, chan_labels};

% PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

% total_dims = sum(window_size > 1) + sum(output_size > 1);
% 
% channel_dims = [2 3];
% 
% indices_cell = make_indices_cell(total_dims, channel_dims, subjects_mat_cell);

SW = PD_sliding_window_pre_post_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin{:});