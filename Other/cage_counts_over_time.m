scriptsDir      = 'C:\Users\alexd\OneDrive - University of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts';
addpath(genpath(scriptsDir));
cd('C:\Users\alexd\Downloads\drive-download-20221213T133713Z-001')
files = dir('down*');
%%

for i = 1 : length (files)
    i
    [~, ~, fExt] = fileparts(files(i).name);
    if strcmp(lower(fExt),'.csv') 
        data = readtable(files(i).name,'ReadVariableNames',1);
    elseif strcmp(lower(fExt),'.xlsx') 
        [num,txt,raw] = xlsread(files(i).name);
        data = table(num,txt(2:end,:));
    end
    % get total cage count
    TCC(i,1)            = size(data,1);
    idx                 = strfind(files(i).name,'-'); idx = idx(1);
    dates{i,1}          = files(i).name(idx-4:idx+5);
    dates_format(i,1)   = datetime(files(i).name(idx-4:idx+5),'InputFormat','yyyy-MM-dd');   
end
%%
% get x tick locations spread according to date
% NumDays = daysact(dates_format(1),dates_format(end))
% months(dates_format(1),dates_format(end))

x = 1;

for i = 2 : size(dates_format,1)
    NumDays = daysact(dates_format(i-1),dates_format(i))
    x(i) = x(i-1) + NumDays
end
%%
hold off
t = table(dates_format,TCC);
t = sortrows(t,'dates_format');
plot(x,t{:,2},'-xk','LineWidth',1,'MarkerSize',10)
xticks(unique(x));
xtickangle(315)
xticklabels(string(unique(t{:,1})))
aesthetics
set(gca,'XMinorTick','off','YMinorTick','off');

MTHS= months(string(dates_format(1)),string(dates_format(end)));
n = 2; % an elbow for every n months
deg = ceil(MTHS/n); % an elbow for every n months
p = polyfit(1:size(t,1), t{:,2}, deg);
v = polyval(p, 1:size(t,1));
hold on
plot(linspace(x(1),x(end),length(x)), v, '--r','LineWidth',5,'Color',.5*[1 1 1])
ylabel('Weekly total NOG cage count')

l = legend('Data',[num2str(deg),' degree polynomial',newline,'(1 elbow per ',num2str(n),' month/s)'],'box','off')
title(["Mouse cage numbers trend from",strcat( string(dates_format(1))," to ",string(dates_format(end)),...
    " (",string(MTHS)," months)" ) ])