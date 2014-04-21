function PD_data_compute_power

close('all')

sampling_freq=2000;
signal_length=sampling_freq*6;

% Parameters for FFT & band power.

f=sampling_freq*[0:signal_length/2]/(signal_length);
f=f(f<=200);
no_freqs=length(f);

spec_format=make_format(no_freqs+1,'f');

band_limits=[.1 4; 4 8; 10 13; 13 20; 20 50; 50 90; 90 120; 125 175];
band_labels={'delta','theta','alpha','low-beta','beta-low-gamma','mid-gamma','high-gamma','HFOs'};

[no_bands,~]=size(band_limits);

band_indices=cell(no_bands,1);
band_freq_labels=cell(no_bands,1);
for i=1:no_bands
    band_indices{i}=find(f>=band_limits(i,1) & f<=band_limits(i,2));
    band_freq_labels{i}=[band_labels{i},num2str(band_limits(i,1)),'-',num2str(band_limits(i,2))];
end

BP_format=make_format(no_bands+1,'f');

% Parameters for blocking out line noise.

line_noise_limits=[59 61; 119 121];

[no_stops,~]=size(line_noise_limits);

stop_indices=cell(no_stops,1);
for i=1:no_stops
    stop_indices{i}=find(f>=line_noise_limits(i,1) & f<=line_noise_limits(i,2));
end

% Getting filenames.

epoch_list='PD_min_epochs.list';
fivemins_list='PD_5min_epochs.list';

[mins,epochs,epoch_names]=textread(epoch_list,'%d%d%s%*[^\n]');
fivemins=textread(fivemins_list,'%d%*[^\n]');

total_epochs=length(epochs);

% Opening files to save collected data.

all_filename=['ALL_',epoch_list(1:end-length('_epochs.list'))];

measures={'spec_pow','band_pow'};
no_measures=length(measures);

fid_pds=fopen([all_filename,'_pds.txt'],'w');

fid_vec=zeros(no_measures,1);
for i=1:no_measures
    fid_vec(i)=fopen([all_filename,'_',measures{i},'.txt'],'w');
end

for i=1:no_bands
    fprintf(fid_vec(2),'%s\t',band_freq_labels{i});
end
fprintf(fid_vec(2),'%s\n','');

spec_all=zeros(total_epochs,no_freqs);
BP_all=zeros(total_epochs,no_bands);

for j=1:total_epochs
    
    fprintf(fid_pds,'%d\t%d\t%d\n',mins(j),fivemins(j),epochs(j));
    
end

%% COMPUTING POWER & BAND POWER BY EPOCH.

parfor j=1:total_epochs
    
    BP=zeros(1,no_bands);
    
    local_band_indices=band_indices;
    
    local_stop_indices=stop_indices;
    
    epoch_name=char(epoch_names(j));
    
    data=load(epoch_name);
    data=detrend(data);
    
    data_hat=pmtm(data);
    
    for i=1:no_stops
        
        data_hat(local_stop_indices{i})=nan;
        
    end
    
    %     data_spec=nanmean(reshape(data_hat(2:end),8,(length(data_hat)-1)/8));
    
    spec_all(j,:)=data_hat(1:no_freqs);
    
    for i=1:no_bands
        
        BP(i)=nansum(data_hat(local_band_indices{i}));
        
    end
    
    BP_all(j,:)=BP;
    
end

%% SAVING POWER & BAND POWER.

fprintf(fid_vec(1),spec_format,[epochs spec_all]');
fprintf(fid_vec(2),BP_format,[epochs BP_all]');

fclose('all');

save([all_filename,'_spec.mat'],'spec_all')
save([all_filename,'_BP.mat'],'band_limits','band_labels','BP_all')