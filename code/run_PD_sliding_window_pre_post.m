%% General stuff

load('seven_bands')

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted')

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

chan_labels = {'Striatum', 'M1'};

filename = 'STR_w_M1';

window_length = 150;

sliding_window_cell{1} = window_length*500*[1 1];

%% Granger causality analysis:

sliding_window_cell{2} = [2 2];

function_handle = @mvgc_analysis; function_name = function_handle;

varargin = {[], '', 1};

% %% PMTM:
% 
% sliding_window_cell{2} = [1 1];
% 
% function_handle = @pmtm;
% 
% varargin = {150, [], 500};

% %% AR spectrum:
% 
% sliding_window_cell{2} = [1 1]; % {[500 500], [1 1]}; % 
% 
% function_handle = @mvgc_analysis; function_name = get_fname(function_handle);
% 
% varargin = {[], '', 3};

% %% AR cross-spectrum:
% 
% sliding_window_cell{2} = [2 2]; % = {[150*500 150*500], [2 2]};
% 
% function_handle = @mvgc_analysis;
% 
% varargin = {[], '', 3};

% %% Running analysis.
% 
% PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

PD_sliding_window_baseline(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

% %% Loading & concatenating analysis.
% 
% [SW, indices_cell] = PD_sliding_window_pre_post_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin{:});
% 
% %% Importing into xPlt.
% 
% SW_xPlt = PD_sliding_window_pre_post_xPlt(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin{:});
% 
% SW_xPlt.getaxisinfo
    
% %% Plotting.
% 
% PD_sliding_window_pre_post_xPlt_plot(function_handle, sliding_window_cell, data_labels_struct, filename, .1, varargin{:})