function indices_cell = make_PD_indices_cell(subjects_mat_cell,function_name,...
    sliding_window_cell, no_windows, output_size, varargin)

no_windows(no_windows == 1) = [];

wdims = length(no_windows);
    
output_size(output_size == 1) = [];

number_odims = length(output_size);

total_dims = 1:(number_odims + wdims);

output_dims = total_dims(1:number_odims);

window_dims = total_dims((number_odims + 1):end);

switch function_name
    
    case 'mvgc_analysis'
        
        switch varargin{end}
            
            case 1
                
                channel_dims = output_dims([2 3]);
                
            case 3
                
                if sliding_window_cell{2}(1) == 1
                    
                    if sliding_window_cell{1}(1) == 150*500
                    
                        channel_dims = window_dims(1);
                    
                    else
                    
                        channel_dims = window_dims(2);
                    
                    end
                        
                elseif sliding_window_cell{2}(1) == 2
                    
                    channel_dims = output_dims([2 3]);
                    
                end
                
        end
        
    case 'pmtm'
        
        if sliding_window_cell{1}(1) == 150
            
            channel_dims = window_dims(1);
            
        else
                    
            channel_dims = window_dims(2);
        
        end
            
end

indices_cell = make_indices_cell(total_dims(end), channel_dims, subjects_mat_cell);

end