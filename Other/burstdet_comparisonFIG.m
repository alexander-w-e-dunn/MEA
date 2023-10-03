
for i=1:length(burstData)
    x(i) = burstData(i).array_burstRate ;
    y(i) = RS_burstData(i).totalBurstRate;
end
[sortedValue_X , X_Ranked] = sort(x,'descend');
[sortedValue_Y , Y_Ranked] = sort(y,'descend');

f1= figure; scatter(X_Ranked,Y_Ranked)
% histogram(x)
% histogram(y)
aesthetics
% rank scores dont correlate which suggests that the burst detection method
% chosen would affect results
xlabel('ISIn method rank')
ylabel('RS method rank')
title('DIV28 MPT cultures array bursting rate (within electrodes)')