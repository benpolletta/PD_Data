%% General stuff

load('seven_bands')

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted')

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

filename = 'STR_w_M1';

sliding_window_cell{1} = 150*[500 500];

%% Granger causality analysis:

sliding_window_cell{2} = [2 2];

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
% varargin = {[], '', 4};

% %% AR cross-spectrum:
% 
% sliding_window_cell = {[150*500 150*500], [2 2]};
% 
% function_handle = @mvgc_analysis;
% 
% varargin = {[], '', 4};

% PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})