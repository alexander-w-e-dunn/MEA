% batch get surrogate connectivity matrices

clear all; close all
%if got spikes already go to spikes folder:
cd 'D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'
% files = dir('*FTD*slice*mSpikes_3.mat');
% files = files(~contains({files.name}, 'TTX','IgnoreCase',true));

% files = dir('*FTD*mSpikes_3.mat');  % where your .mat files are
% files = files(~contains({files.name}, 'ttx'));
% files = files(~contains({files.name}, 'TTX'));
% files = files(~contains({files.name}, 'adjM'));
% files = files(~contains({files.name}, '191210'));
% files = files(~contains({files.name}, 'stim'));
% files = files(~contains({files.name}, 'cleaned'));
% files = files(29);

files = dir('*MPT*cSpikes_L-0.0627_RF1.mat');
files = files(~contains({files.name}, 'TTX','IgnoreCase',true));
files = files(~contains({files.name}, 'DIV07'));
files = files(~contains({files.name}, 'DIV4'));
files = files(~contains({files.name}, 'DIV5'));

sampling_fr = 25000;
% sync_win_s = 0.175;
sync_win_s = 0.05;
method = 'tileCoef';

progressbar('files','iterations','electrodes')
for file = 1:length(files)
    
    tic
    %% get spikes
        data=load(files(file).name,'*Spikes','channels','thresholds'); %may need to at try;catch;end loop if channels var not saved with spikes; could load from mat files
    try
        mSpikes=data.mSpikes;
        disp('loaded mSpikes')
    catch
        cSpikes=data.cSpikes;
        disp('loaded cSpikes')
    end
    
    try
        channels=data.channels;
    catch
        channels=load(strcat(files(file).name(1:18),'.mat'),'channels'); %some spike mats dont have channels var saved
        channels=channels.channels;
        disp('channels loaded from dat file')
    end
    
    try
        thresholds = data.thresholds;
        disp('loaded thresholds')
    catch
    end
    clear data
    
    try
        spikeMat=full(mSpikes);
        clear mSpikes
    catch
        spikeMat=full(cSpikes);
        clear cSpikes
    end
 %% synthesise spike trains and get adjM
    num_iterations = 10;
    for k = 1:num_iterations %for running multiple randomisation iterations
        disp(strcat(num2str(1+num_iterations-k),'_iterations_remaining'))
%         tic
        %%%%%%%%%%%%%%%%%% synthesize trains %%%%%%%%%%%%%%%%%%
        for elec = 1:length(spikeMat(1,:))
            %{
 approach:
                take each 1 ms window and randomly shuffle them
                get a vector of indices for each 1 ms period
                randomly shuffle this
                then resort the spike train using this index
            %}
            
            down_factor = 25;
            fs = 25000;
            fr = sum(spikeMat(:,elec))/(length(spikeMat(:,elec))/fs); %firing rate %divide by seconds
            tSim = length(spikeMat(:,elec))/fs; %duration of simulation in s
            nTrials = 1;%number of trials
            dt = (1/fs)*down_factor; %time bins; set to 1 ms currently
            [rand_spike_vec, tVec] = poissonSpikeGen(fr, tSim, nTrials,dt);
            rand_spikeMat(elec,:) = abs(rand_spike_vec);
            progressbar(file/length(files),k/num_iterations,elec/length(spikeMat(1,:)))
        end
        rand_spikeMat = rand_spikeMat';
        
        %plot to compare spike mat spike counts and rand spike
        %counts over x iterations
        % a(k,:) = sum(spikeMat) - sum(rand_spikeMat);
        
        % %%%%%%%%%%%%%%%%%%% get ctrl adjM    %%%%%%%%%%%%%%%%     
        rec_length_samps = length(rand_spikeMat);
        rec_length_s = rec_length_samps/(sampling_fr/down_factor);
        num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
        ds_rate = length(rand_spikeMat)/length(spikeMat); %scalar for time_window %downsampling rate
        sync_win = sync_win_s * ds_rate;
        adjM2 = getAdjM(rand_spikeMat, method, 0,sync_win); %0 means no downsampling; 0.05 is sync window in s
        
        adjM2 = adjM2 - eye(size(adjM2));
        adjM2(find(isnan(adjM2)))=0;
        %save every random adjM
        adjM2s(:,:,k) = adjM2;
        
        rand_spikeMat_sparse = sparse(rand_spikeMat);
        rand_spikeMat_sparse_all{k,1} = rand_spikeMat_sparse;        
        
        clear adjM2 rand_spikeMat rand_spikeMat_sparse
        toc
        progressbar(file/length(files),k/num_iterations,elec/length(spikeMat(1,:)))
        
    end %for running multiple randomisation iterations
    %                    adjM2 = mean(adjM2s,3);%get mean of random adjMs
%     adjM2 = adjM2s(:,:,1);%get first random adjM
    
    % save control adjMs and spike matrices
    fileNameSpikes = strcat(files(file).name(1:end-4), '_CTRL_spikeMat.mat');
    fileNameadjM = strcat(files(file).name(1:end-4), '_CTRL_adjM_',num2str(sync_win_s),'.mat');
    save(fileNameSpikes, 'rand_spikeMat_sparse_all','-v7.3','channels');
    save(fileNameadjM, 'adjM2s','channels','-v7.3');
    disp('saved first ctrl mat')
    
    clear adjM2 adjM2s rand_spikeMat rand_spikeMat_sparse rand_spikeMat_sparse_all
    progressbar(file/length(files),k/num_iterations,elec/length(spikeMat(1,:)))
    toc
end

