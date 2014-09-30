%run /project/crc-nak/brpp/startup

%cd /project/crc-nak/brpp/PD_Data/

%matlabpool open 8

PD_beta_epochs_GC('st_m1_6OHDA_subjects.mat',7,2,333,5000,2*8*2 + 2)

PD_beta_epochs_GC('st_m1_6OHDA_subjects.mat',7,2,333,20000,0)

PD_beta_epochs_GC('st_m1_subjects.mat',7,2,333,5000,0)

PD_beta_epochs_GC('st_m1_subjects.mat',7,2,333,20000,0)

PD_beta_epochs_GC('st_stn_subjects.mat',7,2,333,5000,0)

PD_beta_epochs_GC('st_stn_subjects.mat',7,2,333,20000,0)

matlabpool close
