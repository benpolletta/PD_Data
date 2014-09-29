%run /project/crc-nak/brpp/startup

cd /project/crc-nak/brpp/PD_Data/

%matlabpool close force local

%matlabpool open local

%PD_beta_epochs_rel_infusion('st_m1_subjects.mat',7,2,1000,5000)

%PD_beta_epochs_rel_infusion('st_stn_subjects.mat',7,2,1000,5000)

%PD_beta_epochs_coh_mtm('st_m1_subjects.mat',7,2,1000,5000,2*1)

%PD_beta_epochs_coh_mtm('st_stn_subjects.mat',7,2,1000,5000,2*1)

%run_PD_coherence('st_m1_subjects.mat',7,2,1000,20000,2*1)

%run_PD_coherence('st_stn_subjects.mat',7,2,1000,20000,2*1)

run_PD_coherence('st_m1_6OHDA_subjects.mat',7,2,1000,20000,2*1)

run_PD_coherence('st_m1_subjects.mat',7,2,2000,20000,2*2)

run_PD_coherence('st_stn_subjects.mat',7,2,2000,20000,2*2)

run_PD_coherence('st_m1_6OHDA_subjects.mat',7,2,2000,20000,2*2)

matlabpool close
