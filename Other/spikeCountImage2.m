% create colour coded image of spike counts 

%% get the spike matrix for that culture
%method = 'Tim'; multiplier = 12; L=0;

%CombineElectrodes;     %run this, see script to change which culture you're analysing
%saves matrix of signals for each electrode
%load the data; may need to change directory
%load allchannels.mat

%spikeMatrix=getSpikeMatrix(data, method, multiplier); %get spikes

spikeCountVector=full(sum(tSpikes)); %note channels are in numerical order, smallest first

    channelMat = [NaN,12:17,NaN;21:28;31:38;41:48;51:58;61:68;71:78;NaN,82:87,NaN]';%note transposed!
    channelVec= [12:17,21:28,31:38,41:48,51:58,61:68,71:78,82:87];
    channelOrder=[47 48 46 45 38 37 28 36 27 17 26 16 35 25 15 14 24 34 13 23 12 22 33 21 ...
        32 31 44 43 41 42 52 51 53 54 61 62 71 63 72 82 73 83 64 74 84 85 75 65 86 76 87 ...
        77 66 78 67 68 55 56 58 57]
    
    for i = 1:length(channelOrder)
        channel_pos(i)=find(channelOrder(i)==channelVec)
    end
    
%determine positions of each channel from 12 to 87
for i=1:length(channelVec)
    channel_mat_pos(i)=find(channelMat==channelVec(i));
end

 
for i =1:60
    spikeCountVec_num_order(channel_pos(i))=spikeCountVector(i);
end

spikeCountMatrix = zeros(8);  
spikeCountMatrix([1,8,57,64])=NaN;
for i = 1:60
    spikeCountMatrix(channel_mat_pos(i))=spikeCountVec_num_order(i)
end 


normalisedSpikeMat = round((1/max(max(spikeCountMatrix)))*spikeCountMatrix,3);
heatmap(normalisedSpikeMat) %v postively skewed so only largest spike counts are visible

logSpikeMat=log(spikeCountMatrix); %note this ln the nat. log. (base=e)
logSpikeMat([1,8,57,64])=NaN;
heatmap(logSpikeMat,'CellLabelColor','none') %plot log spikes
%aesthetics
%cannot add labels with channel numbers to each box!
%simply overlay in powerpoint - see robinson presentation
%insert shape; fill=heatmap from clipboard(use screen snip of matlab figure)
%increase transparency

%define own log:    log9(x) = log(x) / log(9); 
%or could just use log10 command
log100SpikeMat=log(spikeCountMatrix)/log(100);
log100SpikeMat([1,8,57,64])=NaN;
heatmap(log100SpikeMat) %plot log spikes


%% go through each file and get spike matrix and sum them up 


