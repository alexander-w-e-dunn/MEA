%
%
% INPUTS:
%           
%           spikeMatrix:    n samples x n channels, binary
%                           sparse spike matrix runs faster
%                           sparse(spikeMatrix)
%           lag:            synchronicity window (s)
%           fs:             sampling frequency
%           parallel:       if 1, use parallel computing toolbox
%           

function adjM = sttc_fcn(spikeMatrix,lag,fs,parallel)

numChannel = size(spikeMatrix,2);
combChannel = nchoosek(1:numChannel, 2);
A = zeros(1, length(combChannel));
adjM = NaN(numChannel, numChannel);

% make the matrix sparse so that find function runs more quickly
if ~issparse(spikeMatrix)
    spikeMatrix = sparse(spikeMatrix);
end

if parallel == 1
    parfor i = 1:length(combChannel)
        dtv = lag; % [s]
        % spike_times_1 = double(spikeTimes{combChannel(i,1)}.bior1p5/fs);
        % spike_times_2 = double(spikeTimes{combChannel(i,2)}.bior1p5/fs);
        spike_times_1 = find(spikeMatrix(:,combChannel(i,1)) == 1) ./ fs;
        spike_times_2 = find(spikeMatrix(:,combChannel(i,2)) == 1) ./ fs;
        N1v = int64(length(spike_times_1));
        N2v = int64(length(spike_times_2));
        dtv = double(dtv);
        % Time = double([0 spikeDetectionResult.params.duration]);
        Time = [ 0    size(spikeMatrix,1)/fs ];
        tileCoef = sttc_jeremi(N1v, N2v, dtv, Time, spike_times_1, spike_times_2);
        row = combChannel(i,1);
        col = combChannel(i,2);
        A(i) = tileCoef; % Faster to only get upper triangle so might as well store as vector
    end
else
    for i = 1:length(combChannel)
        dtv = lag; % [s]
        % spike_times_1 = double(spikeTimes{combChannel(i,1)}.bior1p5/fs);
        % spike_times_2 = double(spikeTimes{combChannel(i,2)}.bior1p5/fs);
        spike_times_1 = find(spikeMatrix(:,combChannel(i,1)) == 1) ./ fs;
        spike_times_2 = find(spikeMatrix(:,combChannel(i,2)) == 1) ./ fs;
        N1v = int64(length(spike_times_1));
        N2v = int64(length(spike_times_2));
        dtv = double(dtv);
        % Time = double([0 spikeDetectionResult.params.duration]);
        Time = [ 0    size(spikeMatrix,1)/fs ];
        tileCoef = sttc_jeremi(N1v, N2v, dtv, Time, spike_times_1, spike_times_2);
        row = combChannel(i,1);
        col = combChannel(i,2);
        A(i) = tileCoef; % Faster to only get upper triangle so might as well store as vector
    end
end

% Vector -> matrix
for i = 1:length(combChannel)
    row = combChannel(i,1);
    col = combChannel(i,2);
    adjM(row, col) = A(i);
    adjM(col, row) = A(i);
end


end