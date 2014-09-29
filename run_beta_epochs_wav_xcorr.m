run /project/crc-nak/brpp/startup

cd /project/crc-nak/brpp/PD_Data/

matlabpool open local

PD_beta_epochs_wav_xcorr('st_m1_subjects.mat',7,2,2000,20000,8:2:32)

PD_beta_epochs_wav_xcorr('st_stn_subjects.mat',7,2,2000,20000,8:2:32)

matlabpool close