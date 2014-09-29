run /project/crc-nak/brpp/startup

matlabpool open local

% PD_beta_epochs_xcorr('st_m1_subjects.mat',7,2,333,20000,'_A')
% 
% PD_beta_epochs_xcorr('st_stn_subjects.mat',7,2,333,20000,'_A')
% 
% PD_beta_epochs_xcorr('st_m1_subjects.mat',7,2,333,20000,'_P')
% 
% PD_beta_epochs_xcorr('st_stn_subjects.mat',7,2,333,20000,'_P')

PD_beta_epochs_xcorr('st_m1_subjects.mat',7,2,333,20000,'_F')

PD_beta_epochs_xcorr('st_stn_subjects.mat',7,2,333,20000,'_F')

matlabpool close