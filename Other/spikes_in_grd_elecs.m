function percent_correct = spikes_in_grd_elecs(method,parameter)

% script to check no spikes found in ref electrode
    %inputs: the spike detection method ('tim','manuel',or 'cwt', and chosen parameter of threshold
    %or L e.g. 0.2507 if cwt method or 14 if tim method)

%% change to folder with output of checkSpikes.m
if strcmp(method,'cwt') %string compare function compares two inputs, returns logical 1 if same
    if parameter == 0.2507
        cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.2507'
    elseif parameter == 0.0627
        cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\cSpikes.0627'
    else
        display('check parameter folder/file exists')
    end
elseif strcmp(method,'tim')
    if parameter == 14
       cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\tSpikes14'
    elseif parameter == 12
       cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\tSpikes12'
    else 
       display('check parameter folder/file exists')
    end
elseif strcmp(method,'manuel')
    if parameter == 7.5
       cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mSpikes7.5'
    elseif parameter == 5.5
       cd 'D:\MECP2_2019_AD\Data_To_Use\2.4.2.TopCultures\mSpikes5.5'
    else 
       display('check parameter folder/file exists')
    end
else
    display('check method input is correct')
end

%% load output file
datafile = dir('*stats.mat'); 
load(datafile.name,'output');
disp(strcat('analysing:_',datafile.name));

for i = 1:length(output)
    spikes=output(i).spikes; %count of spikes in each elec
    refspikes=spikes(15);
    %this is the ref electrode
    %instead could create loop that checks the filename and then searches a
    %file containing recording names and their corresponding list of grounded
    %electrodes (necessary because these differ between recordings)
    %then set refspikes index to one of the grounded electrodes
    ref_out(i)=refspikes;
end 

badfiles=find(ref_out>0); %files with spikes in ref elec
correct= length(output)-length(badfiles);               %count of files with spikes in ref elec
percent_correct = (correct/length(output))*100;%per cent of the time, the method correctly said 0 in ref channel
disp(strcat('Correct___',int2str(percent_correct),'%',' of the time'));

if percent_correct < 100
    disp(strcat(int2str(length(badfiles)),' ref electrodes have spikes'))
    disp(strcat(int2str(mean(ref_out)),' spikes on mean average across all ref elecs'))
    disp(strcat(int2str(max(ref_out)),' spikes was the maximum detected across all ref elecs'))
else
end

%change dir back to function location to run again to compare methods
cd 'D:\MECP2_2019_AD\Scripts_and_Output';

end