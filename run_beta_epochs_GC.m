%run /projectnb/crc-nak/brpp/startup

%cd /projectnb/crc-nak/brpp/PD_Data/

matlabpool open 8

PD_beta_epochs_GC('st_m1_6OHDA_subjects.mat',7,2,333,5000,8*2 + 3*2)

PD_beta_epochs_GC('st_m1_6OHDA_subjects.mat',7,2,333,20000,0)

PD_beta_epochs_GC('st_m1_subjects.mat',7,2,333,5000,0)

PD_beta_epochs_GC('st_m1_subjects.mat',7,2,333,20000,0)

PD_beta_epochs_GC('st_stn_subjects.mat',7,2,333,5000,0)

PD_beta_epochs_GC('st_stn_subjects.mat',7,2,333,20000,0)

matlabpool close
