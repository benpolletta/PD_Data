%% General stuff

load('seven_bands')

data_labels_struct = init_data_labels(freqs, no_cycles, bands, 'data_field', 'data_subtracted')

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

chan_labels = {'Striatum', 'M1'};

filename = 'STR_w_M1';

if ~exist('window_length', 'var'), window_length = 150; end

sliding_window_cell{1} = window_length*data_labels_struct.sampling_freq{1}*[1 1];

%% Specific stuff.

if ~exist('analysis', 'var'), analysis = 'granger'; end

switch analysis
    
    case 'granger'
        
        % Granger causality analysis:
        
        sliding_window_cell{2} = [2 2];
        
        function_handle = @mvgc_analysis; function_name = function_handle;
        
        varargin = {[], '', 1};
        
        shuffle_struct = init_struct({'shuffle_dims', 'no_shuffles', 'combine_pairs', 'varargin'},...
            {1, 1000, @combine_dimension, {2, [1 2]}});
        
    case 'pmtm'
        
        % PMTM:
        
        sliding_window_cell{2} = [1 1];
        
        function_handle = @pmtm_detrend;
        
        varargin = {window_length, [], data_labels_struct.sampling_freq{1}};
        
    case 'arspec'
        
        % AR spectrum:
        
        sliding_window_cell{2} = [1 1]; % {[500 500], [1 1]}; %
        
        function_handle = @mvgc_analysis; function_name = get_fname(function_handle);
        
        varargin = {[], '', 3};
        
    case 'arxspec'
        
        % AR cross-spectrum:
        
        sliding_window_cell{2} = [2 2]; % = {[150*500 150*500], [2 2]};
        
        function_handle = @mvgc_analysis;
        
        varargin = {[], '', 3};
        
        shuffle_struct = init_struct({'shuffle_dims', 'no_shuffles', 'combine_pairs', 'varargin'},...
            {1, 1000, @combine_dimension, {2, [1 2]}});
        
    case 'pac'
        
        sliding_window_cell{2} = [1 1]; % = {[150*500 150*500], [2 2]};
        
        function_handle = @PAC;
        
        varargin = {data_labels_struct.sampling_freq{1}, [.25:.25:10 11:30], [15:30 32.5:2.5:80 85:5:200]};
        
        shuffle_struct = init_struct({'shuffle_dims', 'no_shuffles', 'combine_pairs', 'varargin'},...
            {1, 1000, @pass_input, {}});
        
end

function_name = get_fname(function_handle);

%% Running analysis.

PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

PD_sliding_window_baseline(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:})

PD_sliding_window_pre_post_shuffle(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, shuffle_struct, varargin{:})

PD_sliding_window_shuffle(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, {'baseline'}, shuffle_struct, varargin{:})

% %% Loading & concatenating analysis.
% 
% [SW, indices_cell] = PD_sliding_window_pre_post_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin{:});

%% Importing into xPlt.

SW_xPlt = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'pre', 'post'}, varargin{:});

SW_xPlt.getaxisinfo
    
SW_xPlt = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'baseline'}, varargin{:});

SW_xPlt.getaxisinfo
    
SW_xPlt = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'pre_shuffles', 'post_shuffles'}, varargin{:});

SW_xPlt.getaxisinfo
    
SW_xPlt = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'baseline_shuffles'}, varargin{:});

SW_xPlt.getaxisinfo

%% Plotting.

whos = {'baseline', 'shuffle'};

hows = {'', 'subtract', 'zscore'};

for w = 1:length(whos)
    
    for h = 1:length(hows)

        norm_struct = struct('who', whos{w}, 'how', hows{h});

        PD_sliding_window_pre_post_xPlt_plot(function_handle, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})

        PD_sliding_window_pre_post_xPlt_bands_plot(function_name, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})
        
    end
    
end