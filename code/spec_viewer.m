datafiles=dir('*.mat');

no_files=length(datafiles);

for f=1:no_files
    
    load(datafiles(f).name)
    
    data=decimate(data,10);
    data=decimate(data,2);
    
    data=detrend(data,'linear');
    
    sampling_rate=1000;
    
    no_secs=floor(length(data)/sampling_rate);
    
    [first_pmtm,f]=pmtm(data(1:sampling_rate),[],[],sampling_rate);
    
    f(f>200)=[];
    
    pmtm_length=length(f);
    
    spectrum=zeros(no_secs,pmtm_length);
    
    spectrum(1,:)=first_pmtm(1:pmtm_length);
    
    for s=1:no_secs
        
        spec_temp=pmtm(data([1:sampling_rate]+(s-1)*sampling_rate),[],[],sampling_rate);
        
        spectrum(s,:)=spec_temp(1:pmtm_length);
        
    end
    
    spec_mean=ones(size(spectrum))*diag(nanmean(spectrum));
    spec_std=ones(size(spectrum))*diag(nanstd(spectrum));
    spec_norm=(spectrum-spec_mean)./spec_std;
    
    t=[1:sampling_rate]/sampling_rate;
    
    figure()
    
    imagesc(t,f,spec_norm)
    
end