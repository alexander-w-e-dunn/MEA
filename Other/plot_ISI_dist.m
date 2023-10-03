%% get ISIs; plot CDF function and theoretical poission CDF function based
% on the mean ISI

% mu = 22;
% sigma = 1;
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% % Define the input vector x to contain the values at which to calculate the cdf.
% 
% x = linspace(22-5,22+5,100);
% % Compute the cdf values for the standard normal distribution at the values in x.
% 
% y = cdf(pd,x);
% figure;
% plot(y)
% 
% %% Create a Poisson distribution object with the rate parameter, ?, equal to 2.
% % 
% % figure; plot(y)
% % figure; plot(pdf(pd,x))
% % 
% figure;
% ISIs = normrnd(10,2.5,100,1);
% histogram(ISIs,'Normalization','probability')
% aesthetics
% hold on
% 
% lambda = mean(ISIs);
% pd = makedist('Poisson','lambda',lambda);
% % Define the input vector x to contain the values at which to calculate the cdf.
% maxISI = 30;
% minISI = 1;
% x = [minISI:1:maxISI];
% % xlim([0 10])
% % Compute the cdf values for the Poisson distribution at the values in x.
% 
% y = cdf(pd,x);
% 
% plot(pdf(pd,x))
% 
% mu = mean(ISIs);
% sigma = std(ISIs);
% pd = makedist('Normal','mu',mu,'sigma',sigma);
% plot(pdf(pd,x))
% % lambda = [1:1:20];
% % y = 1 - poisscdf(2,lambda);
% % pd = makedist('Normal','mu',mu,'sigma',sigma);
% % % Define the input vector x to contain the values at which to calculate the cdf.
% % 
% % x = linspace(22-5,22+5,100);
% % % Compute the cdf values for the standard normal distribution at the values in x.
% % 
% % y = poisscdf(pd,x);
% % figure;
% % plot(y)
% % 

% load data set parameters
load('200127_FTDOrg_GrpD_5B_Slice2_mSpikes_3.mat')
% load('MPT190403_6C_DIV28_cSpikes_L-0.0627_RF1.mat') 
% load('200127_FTDOrg_GrpD_5B_Slice2_mSpikes_3_CTRL_spikeMat.mat')
% mSpikes = rand_spikeMat_sparse_all{1,1};
% if ~exist('fs')
    fs = 25000;
% end
% if ~~exist('rand_spikeMat_sparse_all')
%     fs = 1000;
% end
try
    mSpikes_real = mSpikes;
catch
    mSpikes = cSpikes;
    mSpikes_real = mSpikes;
end
% clear mSpikes
spiketimes = find(full(mSpikes_real(:,1))==1);
ISIs = ( spiketimes(2:end,1) - spiketimes(1:end-1,1) )  / (fs/1000);
% ISIs = log(ISIs);
% spiketimes2 = find(full(mSpikes_synt(:,1))==1);
% ISIs2 = ( spiketimes2(2:end,1) - spiketimes2(1:end-1,1) )  / (1000/1000);

f1 = figure;
h = histogram(sqrt(ISIs),'Normalization','probability')
h = histogram(log10(ISIs),'Normalization','probability')
h = histogram(ISIs,'Normalization','probability','FaceColor','k');  hold on
a = yticks; alen = length(a);
a1 = yticklabels;
yl = ylim;
h2 = histfit(ISIs,100,'gamma')  organoid
% https://uk.mathworks.com/help/stats/histfit.html 
% h2 = histfit(ISIs,5000,'gamma') % Mecp2
% h2(1).Visible = 'off'
ylim(yl);
h2(2).Color = 'r';
h2(2).LineWidth = 3;
h2(2).LineStyle = '--';
h2(1).FaceColor = 0*[1 1 1];
h2(1).EdgeColor = 1*[1 1 1];
b = yticks;
yticks( linspace(min(b), max(b),alen))
yticklabels(a1);

%options: inversegaussian exponential kernel gev

% peak = 1 + h.BinWidth * x(find(h.Values == max(h.Values)));

aesthetics
% xticklabels(xticks.^2)
% xticklabels(10.^xticks)
xticklabels(xticks)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('p');
% ylabel('{p       }')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'FontName','Arial','FontSize',14)
l1 = legend('real','gamma fit','box','off');
f1.Position = [680   728   499   250];
rxlims = xlim;
rylims = ylim;

spikeMat = full(mSpikes);
spikeTimes = findSpikeTimes(spikeMat(:,1), 'seconds', fs);
spikeISI = findISI(spikeTimes);
regularity = getReg(spikeISI, 'gamma', 8000) % for organoid
regularity = getReg(spikeISI, 'gamma', 500) %for mecp2
% ax = gca;
% ax.FontSize = 12
% histogram(ISIs)
clear mSpikes
%% get same for rand mat
load('200127_FTDOrg_GrpD_5B_Slice2_mSpikes_3_CTRL_spikeMat.mat')
mSpikes = rand_spikeMat_sparse_all{1,1};
if ~exist('fs')
    fs = 25000;
end
if ~~exist('rand_spikeMat_sparse_all')
    fs = 1000;
end
mSpikes_real = mSpikes;
spiketimes = find(full(mSpikes_real(:,1))==1);
ISIs = ( spiketimes(2:end,1) - spiketimes(1:end-1,1) )  / (fs/1000);
% ISIs = log(ISIs);
% spiketimes2 = find(full(mSpikes_synt(:,1))==1);
% ISIs2 = ( spiketimes2(2:end,1) - spiketimes2(1:end-1,1) )  / (1000/1000);

f1 = figure;
% h = histogram(sqrt(ISIs),'Normalization','probability')
% h = histogram(log10(ISIs),'Normalization','probability')
% h = histogram(ISIs,'Normalization','probability','FaceColor','k'); hold on
% a = yticks; alen = length(a);
% a1 = yticklabels;
h2 = histfit(ISIs,100,'gamma') % https://uk.mathworks.com/help/stats/histfit.html 
h2(2).Color = 'r';
h2(2).LineWidth = 3;
h2(2).LineStyle = '--';
h2(1).FaceColor = 0.5*[1 1 1];
h2(1).EdgeColor = 1*[1 1 1];
b = yticks;
yticks( linspace(min(b), max(b),alen))
yticklabels(a1);

%options: inversegaussian exponential kernel gev

% peak = 1 + h.BinWidth * x(find(h.Values == max(h.Values)));

aesthetics
% xticklabels(xticks.^2)
% xticklabels(10.^xticks)
xticklabels(xticks)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('p');
% ylabel('{p       }')
% set(get(gca,'ylabel'),'rotation',0)
set(gca,'FontName','Arial','FontSize',14)
l1 = legend('synth.','gamma fit','box','off');
f1.Position = [680   728   499   250];
xlim(rxlims);ylim(rylims);

spikeMat = full(mSpikes);
spikeTimes = findSpikeTimes(spikeMat(:,1), 'seconds', fs);
spikeISI = findISI(spikeTimes);
regularity = getReg(spikeISI, 'gamma', 8000)

%%
X = ISIs;
% X2 = ISIs2;
[f,x_values] = ecdf(X);
% [f2,x_values2] = ecdf(X2);
figure;
J = plot(x_values,f,'k');
hold on;
% J2 = plot(x_values2,f2,'c--');
K = plot(x_values,normcdf(x_values,mean(x_values),std(x_values)),'r--');
lambda = mean(X);
% y = cdf('Poisson',X,lambda);

% y = poisscdf(x_values2,lambda);
% plot(1:length(y),sort(y),'m--');

mu = mean(X);
y = expcdf(x_values,mu);
I = plot(x_values,sort(y),'m--');

% y = sort(y); x = 1:length(y); I = plot(x_values,y(1:length(x_values)),'m--');
% pd = makedist('Poisson','lambda',lambda);
% PoissX = cdf(pd,X);
% I = plot(X,PoissX,'m--');
% PoissX = poisscdf(X,lambda);
% I = plot(x_values,poisscdf(x_values,lambda),'m--');
set(J,'LineWidth',2);
% set(J2,'LineWidth',2);
set(K,'LineWidth',2);
set(I,'LineWidth',2);
legend([J K I],'Empirical CDF','Standard Normal CDF','Exponential CDF','Location','SE');
aesthetics
% xticklabels(xticks.^2)
% xticklabels(10.^xticks)
xticklabels(xticks)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('{p       }')
set(get(gca,'ylabel'),'rotation',0);
set(gca,'FontName','Arial','FontSize',14);

%% 
mean(sqrt(ISIs))

lambda = mean(sqrt(ISIs));
% lambda = peak;
pd = makedist('Poisson','lambda',lambda);
% Define the input vector x to contain the values at which to calculate the cdf.
% maxISI = max(sqrt(ISIs));
% minISI = min(sqrt(ISIs));
x = [1:1:max(xticks)];
% xlim([0 10])
% Compute the cdf values for the Poisson distribution at the values in x.

%y = cdf(pd,x);
hold on

mu =  mean(sqrt(ISIs));
% mu = peak;
 sigma = std(sqrt(ISIs));

pdn = makedist('Normal','mu',mu,'sigma',sigma);
plot(pdf(pdn,x)/2,'LineWidth',2,'Color','b','LineStyle',':')

plot(pdf(pd,x)/2,'LineWidth',2,'Color','r','LineStyle',':')


%cdf plot 

figure
plot(cdf(pdn,x),'LineWidth',2,'Color','b','LineStyle',':')
hold on
plot(cdf(pd,x),'LineWidth',2,'Color','r','LineStyle',':')
aesthetics
xticklabels(xticks.^2)
xtickformat(    )
xlabel('Inter-spike interval (ms)')
ylabel('Cumulative p')
set(get(gca,'ylabel'),'rotation',90)
set(gca,'FontName','Arial','FontSize',14)

c =  cdfplot(sqrt(ISIs));  
grid off
c.LineWidth = 2;
c.Color = 'k';
c.Parent.XLabel.String = 'Inter-spike interval (ms)';
legend('Poisson CDF','Normal CDF','Empirical CDF','Location','best');
% sq_sorted = sort(sqrt(ISIs));
% cdf(sq_sorted)
