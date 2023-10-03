%% plot filtered trace of spikes overlaid and peaks aligned
function spike_overlaysingle_fcn(spikeTrain,finalData,method,parameter,refPeriod_ms)

multiplier          = parameter;
L                   = parameter;


%find channel index (note this is the column number of the desired data,
%not the MCS XY coordinate)
%not MCS ID but which column the desired E is in
electrode_to_plot   = 1;
%correction if it's a draw (i.e. >1 electrodes has max no. spikes)
if  length(electrode_to_plot)>1
    electrode_to_plot = electrode_to_plot(1);
end


    
    %% get spike times and plot trace
    
    sp_times=find(spikeTrain==1);
    n_spikes_to_plot=50;
    %added correction if num spikes is fewer than desired number:
    if  sum(spikeTrain) < n_spikes_to_plot
        n_spikes_to_plot = sum(spikeTrain);
    end
    
    for i=1:n_spikes_to_plot
        sp_peak_time=find(finalData(sp_times(i):sp_times(i)+25)==min(finalData(sp_times(i):sp_times(i)+25)));%note changed min to max as with filtered data spike is +ve
        %plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25));
        
        % added correction for if there is an early, such that spike time -
        % minus 25 samples is negative, it will not plot
        % or if any part of the trace to plot is above 50 or below -50, it
        % will exclude this trace from the plot and the average
        % calculation
        if (sp_times(i))+sp_peak_time-25 < 1 | ...
                ~isempty(find(...
                finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25) > 50)) | ...
                ~isempty(find(...
                finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25) < -50))
            % don't plot; don't add average, first line of array will be
            % excluded by MATLAB automatically
        else
            p1 = plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
                'Color',0.4*[1 1 1],'LineWidth',0.5); %all grey
            p1.Color(4) = 0.2;
%             plot(finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25),...
%                 'LineWidth',0.5);
            hold on
            all_trace(:,i)=finalData((sp_times(i))+sp_peak_time-25:(sp_times(i))+sp_peak_time+25);
        end
    end
    
    ave_trace=(sum(all_trace'))/(length(all_trace(1,:)));
    %plot(ave_trace,'LineWidth',8,'Color','w');
    %plot(ave_trace,'LineWidth',3,'Color','r');
    plot(ave_trace,'LineWidth',1.5,'Color','k'); %black line instead
    hold off
    aesthetics
    
    %change axes to voltage normalised
    %change x axis to time rather than samples
    xticks(linspace(0,50,3));
    xticklabels([linspace(0,50,3)/25]-1);
    xlabel('time relative to negative peak (ms)');
    
    %calibrate y axis
    nearestValue = 50; % i.e. nearest mutliple of 50
    ymax = ceil(max(max(all_trace))/nearestValue)*nearestValue;
    ymin = ceil(min(min(all_trace))/nearestValue)*nearestValue;
    if  abs(ymax) >= abs(ymin)
        yval = abs(ymax);
    else
        yval = abs(ymin);
    end
    
    ylim([-yval yval])
    
    if      yval >= 300
        increment = 100;
    elseif  yval >= 100
        increment = 50;
    else
        increment = 25;
    end
    
    yticks(linspace(-yval,yval,((2*yval)/increment)+1));
    ylabel('filtered signal (\muV)');
    set(gca,'fontsize',16)
    
    % plot the threshold if using threshold-based method
    if      strcmp(method,'Manuel') 
        threshold = mean(finalData(:,1)) - multiplier * std((finalData(:,1)));  
    elseif strcmp(method,'abs')
        threshold = parameter;
    end
    
%     if      strcmp(method,'Manuel') | strcmp(method,'abs')
%         threshvec = ones(length(ave_trace))*threshold;
%         hold on
%         plot(1:length(ave_trace),threshvec,'LineStyle',':','Color','r','LineWidth',2)
        axis off; xlim([-20 71]);ylim([-50 50]);
        hold off
%         %calculate the average and overlay that onto the figure as a red line
%         %surrounded by white
%     else
%     end
    




end