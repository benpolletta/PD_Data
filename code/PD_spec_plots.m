function PD_spec_plots

sampling_freq=2000;
signal_length=sampling_freq*6;

freqs=sampling_freq*[0:signal_length]/signal_length;
freqs=freqs(freqs<=200);

bands=[130 150; 35 90; 12 30];
band_names={'HFO','gamma','beta'};

stops=[58 62; 118 122; 179 181];

[mins,fivemins]=textread('ALL_PD_min_pds.txt','%s%s%*[^\n]');
load('ALL_PD_min_spec.mat');
for i=1:length(mins), dummy{i}=''; end

%% Plots by Minute.

no_pre=5; no_post=14; no_mins=no_pre+no_post+1;
min_labels=cell(no_mins,1); short_min_labels=cell(no_mins,1); min_corder=zeros(no_mins,3);

p=1;
for i=no_pre:-1:1
    min_labels{p}=['Minute ',num2str(i),' Preinjection'];
    short_min_labels{p}=num2str(i);
    min_corder(p,:)=(p-1)*[1 1 1]/(2*no_pre);
    p=p+1;
end
for i=0:no_post
    min_labels{p}=['Hour ',num2str(i),' Postinjection'];
    short_min_labels{p}=num2str(i);
    min_corder(p,:)=i*[0 1 1]/no_post+(no_post-i)*[1 0 1]/no_post;
    p=p+1;
end

cplot_collected_spec_by_categories('Spectral Power','PD_min_spec',freqs,bands,band_names,stops,min_corder,{{''},{''}},{short_min_labels, min_labels},dummy,mins,zscore(spec_all))

%% Plots by 5 Minute Periods.

no_pre=1; no_post=2; no_5mins=no_pre+no_post+1;
fivemin_labels=cell(no_5mins,1); short_5min_labels=cell(no_5mins,1); fivemin_corder=zeros(no_5mins,3);

p=1;
for i=1:no_pre
    fivemin_labels{p}=['Minutes ',num2str(5*i),' to ',num2str(5*(i-1)+1),' Preinjection'];
    short_5min_labels{p}=num2str(-i);
    fivemin_corder(p,:)=(p-1)*[1 1 1]/(2*no_pre);
    p=p+1;
end
for i=0:no_post
    fivemin_labels{p}=['Minutes ',num2str(5*i),' to ',num2str(5*(i+1)-1),' Postinjection'];
    short_5min_labels{p}=num2str(i);
    fivemin_corder(p,:)=i*[0 1 1]/no_post+(no_post-i)*[1 0 1]/no_post;
    p=p+1;
end

cplot_collected_spec_by_categories('Spectral Power','PD_5min_spec',freqs,bands,band_names,stops,fivemin_corder,{{''},{''}},{short_5min_labels, fivemin_labels},dummy,fivemins,zscore(spec_all))

end

