function run_beta_epochs_GC(group_names, win_lengths, smooth_lengths, start_indices)

if isempty(start_indices)
    
    start_indices = zeros(length(group_names)*length(win_lengths)*length(smooth_lengths));
    
elseif start_indices ~= length(group_names)*length(win_lengths)*length(smooth_lengths);
    
    display('The length of optional argument start_indices should be the product of the lengths of group_names, win_lengths, and smooth_lengths.')
    
    return
    
end

run /projectnb/crc-nak/brpp/startup

cd /projectnb/crc-nak/brpp/PD_Data/

s = matlabpool('size');

if s == 0
    
    matlabpool open 8
    
end

index = 1;

for g = 1:group_names
    
    for w = 1:win_lengths
        
        for s = 1:smooth_lengths

            PD_beta_epochs_GC(group_names{g}, 7, 2, win_lengths(w), smooth_lengths(s), start_indices(index))
            
            PD_beta_epochs_GC_collect(group_names{g}, 7, 2, win_lengths(w), smooth_lengths(s))
            
            PD_beta_epochs_GC_plot_group(group_names{g}, 7, 2, win_lengths(w), smooth_lengths(s))
            
        end
        
    end
    
end

if  s == 0
    
    matlabpool close
    
end


