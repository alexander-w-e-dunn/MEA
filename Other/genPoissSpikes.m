 active_spike_mat = spikeMat(:,active_chanIndex);
    for tarly = 1:length(active_spike_mat(1,:))
        %% generate poisson spike train of same rate
        fr = (sum(active_spike_mat(:,tarly)))     /    (length(active_spike_mat(:,tarly))/sampling_fr);%fir rate in Hz
        dt = 1/25000; %time duration of one sample in seconds
        nBins=18000000; %would have to loop through and get firing rate in first e.g. 250,000bins(10s) then next etc
        rsvec=rand(1, nBins) < fr*dt;
        rsvec=double(rsvec)';
        sum(rsvec)
        sum(active_spike_mat(:,tarly))
        %produces similar spike count
        
        %% randomised spike train
        rand_spike_vec=active_spike_mat(:,tarly);
        rand_spike_vec=rand_spike_vec(randperm(length(rand_spike_vec)));
        
        %check distributions
        sts=find(active_spike_mat(:,tarly)==1);
        sts2=sts(2:length(sts));
        ISIsss=sts2-sts(1:end-1);
        ISIsss=ISIsss/25000;
        figure
        hist(ISIsss)
        m=max(ISIsss);
        
        box off
xlabel('ISI (s)')
ylabel('Frequency')
aesthetics
xlim([0 m])
ylim([0 sum(active_spike_mat(:,tarly))])
        
                sts=find(rsvec==1);
        sts2=sts(2:length(sts));
        ISIsss=sts2-sts(1:end-1);
        ISIsss=ISIsss/25000;
        figure
        hist(ISIsss)
        xlim([0 m])
ylim([0 sum(active_spike_mat(:,tarly))])
        
                        sts=find(rand_spike_vec==1);
        sts2=sts(2:length(sts));
        ISIsss=sts2-sts(1:end-1);
        ISIsss=ISIsss/25000;
        figure
        hist(ISIsss)
        xlim([0 m])
ylim([0 sum(active_spike_mat(:,tarly))])
    end 
%cf https://praneethnamburi.com/2015/02/05/simulating-neural-spike-trains/