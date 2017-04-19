function [SW_xPlt, SW_Baseline, SW_Shuffles] = run_PD_sliding_window_pre_post_function(band_index, analysis, window_length, analysis_flag, load_flag, plot_flag)

if nargin == 0, band_index = []; end
if nargin < 2, analysis = []; end
if nargin < 3, window_length = []; end
if nargin < 4, analysis_flag = []; end
if nargin < 5, load_flag = []; end
if nargin < 6, plot_flag = []; end

if isempty(analysis_flag), analysis_flag = [1 1 1 0]; end
if isempty(load_flag), load_flag = 1; end
if isempty(plot_flag), plot_flag = [1 0]; end

%% General stuff

seven_bands = load('seven_bands');

if isempty(band_index), band_index = 4; end
data_labels_struct = init_data_labels(seven_bands.freqs, seven_bands.no_cycles, seven_bands.bands, 'data_field', 'data_subtracted', 'band_index', band_index);

subjects_mat_cell = {'st_m1_subjects.mat', 'st_m1_ali_subjects.mat', 'st_m1_ali2_subjects.mat'};

chan_labels = {'Striatum', 'M1'};

filename = 'STR_w_M1';

if isempty(window_length), window_length = 10; end

sliding_window_cell{1} = window_length*data_labels_struct.sampling_freq{1}*[1 1];

%% Specific stuff.

if isempty(analysis), analysis = 'granger'; end

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

if analysis_flag(1)
    
    PD_sliding_window_pre_post(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:}), end

if analysis_flag(2)
    
    PD_sliding_window_baseline(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, varargin{:}), end

if analysis_flag(3)
    
    PD_sliding_window_pre_post_shuffle(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, shuffle_struct, varargin{:}), end

if analysis_flag(4)
    
    PD_sliding_window_shuffle(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, {'baseline'}, shuffle_struct, varargin{:}), end

% %% Loading & concatenating analysis.
% 
% [SW, indices_cell] = PD_sliding_window_pre_post_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, varargin{:});

%% Importing into xPlt.

if load_flag
    
    SW_xPlt = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'pre', 'post'}, varargin{:});
    
    SW_xPlt.getaxisinfo
    
    SW_Baseline = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'baseline'}, varargin{:});
    
    SW_Baseline.getaxisinfo
    
    SW_Shuffles = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'pre_shuffles', 'post_shuffles'}, varargin{:});
    
    SW_Shuffles.getaxisinfo
    
    SW_BaselineShuffles = PD_sliding_window_load(function_handle, sliding_window_cell, subjects_mat_cell, data_labels_struct, filename, {'baseline_shuffles'}, varargin{:});
    
    SW_BaselineShuffles.getaxisinfo
    
end

%% Plotting.

if plot_flag(1)
    
    whos = {'baseline', 'shuffle'};
    
    hows = {'', 'subtract', 'zscore'};
    
    for w = 1:length(whos)
        
        for h = 1:length(hows)
            
            norm_struct = struct('who', whos{w}, 'how', hows{h});
            
            PD_sliding_window_pre_post_xPlt_plot(function_handle, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})
            
            PD_sliding_window_pre_post_xPlt_bands_plot(function_name, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})
            
            for h1 = 1:length(hows)
                
                norm_struct = struct('who', whos{w}, 'how', {hows{h}, hows{h1}});
                
            end
            
        end
        
    end
    
end

if plot_flag(2)
    
    whos = {'shuffle_baseline', 'baseline_shuffle'};
    
    hows = {'', 'subtract', 'zscore'};
    
    for w = 1:length(whos)
    
        for h = 1:length(hows)
    
            for h1 = 1:length(hows)
    
                norm_struct = struct('who', whos{w});
                
                norm_struct.how = {hows{h}, hows{h1}};
    
                PD_sliding_window_pre_post_xPlt_plot(function_handle, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})
    
                PD_sliding_window_pre_post_xPlt_bands_plot(function_name, sliding_window_cell, data_labels_struct, filename, .1, norm_struct, varargin{:})
    
            end
    
        end
    
    end
    
end

