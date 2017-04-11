% This file takes a file input and generates peaks for it
function pipeline_single(fi_no, dim, min_prominences, clusters, normalized, invert)

    load('STR_M1_subjects.mat')
    
    cd(folders{fi_no})
    
    load([prefixes{fi_no} '_all_channel_data_dec.mat']);

    basetime = basetimes(fi_no);
    
    prefix = prefixes{fi_no};

    PD_data = PD_dec;
    
    fq = sampling_freq;

    [a, b] = size(PD_data);

    if b > a
       PD_data = PD_data';
    end
    
    if invert
        PD_data(:,2) = -PD_data(:,2);
    end
    
    min_secs_apart = .5;
    
    Peak_data = classify_peaks(PD_data(:,dim), fq, min_prominences, min_secs_apart, basetime, clusters, normalized, [prefix '_ch' num2str(dim)]);
    disp(['Done with ' prefix ' channel ' num2str(dim)])
    
    cd ..

end

% best parameters list
% no artifact:
%  mice: 130703, 130709, 130716, 130718 (maybe), 130725, 130822,
%  130830,mice1
% 
% artifact containing:
%  13813 (both channels), 13815 (both channels, but not stereotyped very
%  well), mice2, mice3
% Below is for mouse 13813:
% pipeline_single(6, 100, 3,1,1)
% Below is for mouse 2
% pipeline_single(11, 100, 3,0,0)