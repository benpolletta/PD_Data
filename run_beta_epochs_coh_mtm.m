function run_beta_epochs_coh_mtm(group_names, win_lengths, smooth_lengths)

run /projectnb/crc-nak/brpp/startup

cd /projectnb/crc-nak/brpp/PD_Data/

% s = matlabpool('size');
% 
% if s == 0
%     
%     matlabpool open 8
%     
% end

for g = 1:length(group_names)
    
    for w = 1:length(win_lengths)
        
        for s = 1:length(smooth_lengths)
            
            run_PD_coherence(group_names{g}, 7, 2, win_lengths(w), smooth_lengths(s), 2*win_lengths(w)/1000)

        end
        
    end
    
end

% if  s == 0
%     
%     matlabpool close
%     
% end