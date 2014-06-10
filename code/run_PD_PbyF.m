function run_PD_PbyF

present_dir = pwd;

% load('initial_subjects.mat')

load('st_m1_subjects.mat')

% for fo = 5:length(folders)
%     
%     folder = folders{fo};
%     
%     prefix = prefixes{fo};
%    
%     cd (folder)
%     
%     PD_concatenate_channels(prefix,1)
%     
%     PD_decimate_channels(prefix,1)
% 
%     cd (present_dir)
%     
% end
% 
% PD_bandpass_channels

% PD_beta_epochs(2,round(1000/3),100)

PD_beta_epochs_roseplot_by_datapoint

% PD_beta_epochs_rel_infusion(2,round(1000/3))

PD_beta_epochs_rel_infusion_roseplot_by_datapoint