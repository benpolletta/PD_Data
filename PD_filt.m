load('PD_dec.mat')

slo = 0.05;
plo = 0.1;
phi = 500;
shi = 550;

[n,Wn]=buttord(2*[plo phi]/sampling_freq,2*[slo shi]/sampling_freq,1,70);

[z,p,k]=butter(n,Wn); [sos,g]=zp2sos(z,p,k); h=dfilt.df2sos(sos,g);
[H,W]=freqz(h,nyquist_index);

PD_filtered=filtfilthd(h,PD_dec,'reflect');

saveas('PD_filt.mat','PD_filtered','sampling_freq')