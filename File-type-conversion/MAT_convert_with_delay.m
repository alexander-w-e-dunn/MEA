hours = 2; % number of hours to delay onset of raw to mat conversion
cd 'C:\Users\Windows\Documents\MATLAB'
addpath(genpath('C:\Users\Windows\Documents\MATLAB'))
%input #2 to below timer function is the delay in s
T = timer('StartDelay',hours*60*60,'TimerFcn',@(src,evt)MEAbatchConvert_alex); 
start(T)
