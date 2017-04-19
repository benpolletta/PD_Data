function SW_norm = PD_xPlt_prep_for_norm(SW_norm_struct, SW_xPlt, function_name, pd_names, collapse_struct)
    
    SW_norm = SW_norm_struct.SW_xPlt;
    
    axes_info_struct = SW_norm_struct.axes_info_struct;

    % Taking absolute value if not PAC.
    if ~strcmp(function_name, 'PAC')
        
        SW_norm = xp_abs(SW_norm);
        
    end
        
    % Handling period axis replacement.
    norm_period_dim = SW_norm.findaxis('Period');
    data_period_dim = SW_xPlt.findaxis('Period');
    
    if any(cellfun(@(x) ~isempty(x), regexp(pd_names, 'shuffle', 'once')))
        
        if ~isempty(data_period_dim) && ~isempty(norm_period_dim)
            
            SW_norm.axis(norm_period_dim).values = SW_xPlt.axis(data_period_dim).values;
            
        end
        
    end
    
    if isempty(norm_period_dim) && ~isempty(data_period_dim)
        
        SW_norm = SW_norm.repmat(SW_xPlt.axis(data_period_dim).values, 'Period', data_period_dim);
        
    end
    
    if isempty(data_period_dim) && ~isempty(norm_period_dim)
        
        SW_norm.axis(norm_period_dim).values = pd_names;
        SW_xPlt = SW_xPlt.repmat(pd_names, 'Period', norm_period_dim);
        
    end
    
    % If these are shuffles, then rewriting/unpacking Shuffles axis as Window_Dim_1 axis.
    if any(cellfun(@(x) ~isempty(x), regexp(pd_names, 'shuffle', 'once')))
        
        if ~isempty(SW_norm.findaxis('Shuffles')) % If shuffles are unpacked.
            
            SW_norm.axis(findaxis('Shuffles')).name = 'Window_Dim_1'; % = mean_over_axis(SW_norm, 'Shuffles', collapse_struct.function, collapse_struct.varargin);
            
        else % If shuffles are packed.
            
            SW_norm = SW_norm.unpackDim(axes_info_struct.odims + 1, [], 'Window_Dim_1');
            
        end
        
    end
    
    % If there is a sliding window, taking mean over sliding window for
    % SW_norm values, and replicating over sliding window of SW_xPlt values.
    if ~isempty(SW_norm.findaxis('Window_Dim_1')) % If norm values are unpacked.
        
        SW_norm = mean_over_axis(SW_norm, 'Window_Dim_1', collapse_struct.function, collapse_struct.varargin);
        
    elseif isfield(SW_norm.meta, ['matrix_dim_', num2str(axes_info_struct.odims + 1)])
        
        varargin = collapse_struct.varargin;
        
        varargin{end} = axes_info_struct.odims + 1;
        
        SW_norm.data = cellfun(@(x) eval(collapse_struct.function, x, varargin{:}), SW_norm.data, 'UniformOutput', false);
        
        SW_norm.meta = rmfield(SW_norm.meta, ['matrix_dim_', num2str(axes_info_struct.odims + 1)]);
        
    end
    
    if ~isempty(SW_xPlt.findaxis('Window_Dim_1')) % Unpacking mean value over window dimension of xPlt.
        
        SW_norm = SW_norm.repmat(SW_xPlt.axis(SW_xPlt.findaxis('Window_Dim_1')).values, 'Window_Dim_1', SW_xPlt.findaxis('Window_Dim_1'));
        
    end
        
    SW_norm = SW_norm.alignAxes(SW_xPlt);

end