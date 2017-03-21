function indices_cell = make_indices_cell(total_dims, channel_dims, subjects_mat_cell)

indices_cell = cell(length(subjects_mat_cell), 2);

for s = 1:length(subjects_mat_cell)
    
    subjects_mat_struct = load(subjects_mat_cell{s});
   
    index_cell = cell(1, total_dims);
    
    index_cell(:) = {':'};
    
    [input_index_cell, output_index_cell] = deal(index_cell);
    
    input_index_cell(channel_dims) = {[1 2]};
    
    chan_order = [1 2];
    
    if strcmp(subjects_mat_struct.chan_labels{2}, 'Striatum')
        
        chan_order = fliplr(chan_order);
        
    end
    
    output_index_cell(channel_dims) = {chan_order};
    
    indices_cell(s, :) = {output_index_cell, input_index_cell};
    
end

end