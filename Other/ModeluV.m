% code for modelling and plotting membrane potential fot the purposes of
% MEA data
% Written by Alex Dunn, University of Cambridge, 2020
% This script uses code from www.izhikevich.com by Eugene M. Izhikevich, February 25, 2003
% E. M. Izhikevich, “Simple model of spiking neurons,” IEEE Transactions on Neural Networks, vol. 14, no. 6. pp. 1569–1572, Nov-2003.
clear all
close all
addpath(genpath('D:\MECP2_2019_AD\Scripts_and_Output\S1.2.File_Conversion_Output'))
MEAgraphics(0); MEAgraphics(1)
%% Izhikevich model:
% This section was created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons    Inhibitory neurons
Ne=800;                 Ni=200;
re=rand(Ne,1);          ri=rand(Ni,1);
a=[0.02*ones(Ne,1);     0.02+0.08*ri];
b=[0.2*ones(Ne,1);      0.25-0.05*ri];
c=[-65+15*re.^2;        -65*ones(Ni,1)];
d=[8-6*re.^2;           2*ones(Ni,1)];
S=[0.5*rand(Ne+Ni,Ne),  -rand(Ne+Ni,Ni)];

v=-65*ones(Ne+Ni,1);    % Initial values of v
u=b.*v;                 % Initial values of u
firings=[];             % spike timings

for t=1:1000            % simulation of 1000 ms
  I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
  fired=find(v>=30);    % indices of spikes
  firings=[firings; t+0*fired,fired];
  v(fired)=c(fired);
  u(fired)=u(fired)+d(fired);
  I=I+sum(S(:,fired),2);
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % step 0.5 ms
  v=v+0.5*(0.04*v.^2+5*v+140-u+I); % for numerical
  u=u+a.*(b.*v-u);                 % stability
end;

%% 

figure
subplot(4,1,1:3)
plot(firings(:,1),firings(:,2),'.k');
aesthetics
hold on
plot([1 length(v)],[Ne Ne],'-k','LineWidth',1)
a=gca;
a.XAxis.TickLabels = [];
a.XAxis.Label = [];
a.YLabel.String = 'Neuron number';
a.YTickLabelRotation = 90;
a.FontSize = 12;
yticks([1 Ne length(v)])
yts = yticks;
% text(length(v)+20,length(v),'Inhibitory','Rotation',270,'FontSize',12,'FontName','Arial')
text(length(v)+20,mean([length(v) Ne]),'Inh.','Rotation',0,'FontSize',12,'FontName','Arial')
text(length(v)+20,mean([0 Ne]),'Exc.','Rotation',0,'FontSize',12,'FontName','Arial')
arpos = a.Position;
% an = annotation('textarrow',[arpos(1)-0.03,arpos(1)-0.03],[0.35 0.75],...
%     'String','Neuron number','FontSize',12,'Linewidth',2,...
%     'TextRotation',90);
annotation('textarrow',[arpos(1)-0.03,arpos(1)-0.03],[0.35 0.75],...
    'FontSize',12,'Linewidth',2)
hold off
subplot(4,1,4)
plot(v,'k-','LineWidth',1)
ylim([-90 +40]);
aesthetics
b=gca;
b.YLabel.String = '\muV';
b.FontSize = 12;
b.TickDir = 'out';
a.TickDir = 'out';
b.YTickLabelRotation = 0;
yticks([-65 0 40]);
yticklabels([-65 0 40]);
xlabel('Time (ms)')
f=gcf;
f.Position = [700   281   596   504];

%% 
spikemat = zeros(sum([Ne Ni]));
for i = 1 : length(firings)
    spikemat(firings(i,1),firings(i,2)) = 1;
end
% correlate using sttc and 50 ms sync window
am=getAdjM(spikemat, 'tileCoef', 0, 50/25000); % sync window is in s/25000
% because getAdjM function assumes sample rate is 25000 but sttc function
% requires sync window in s
%%
figure;
subplot(1,2,1);
imagesc(S(length(S):-1:1,:));
cb= colorbar
cb.Location = 'Southoutside';
cb.Label.String = 'Synaptic weight';
ylabel('Neuron (post)');
xlabel('Neuron (pre)');
yts = yticks;
yticks([1 200:200:800])
yticklabels(flip(yts));
text(mean([length(v) Ne])-50,length(v)-length(v)-50,'Inh.','Rotation',0,'FontSize',12,'FontName','Arial')
text(mean([0 Ne]),length(v)-length(v)-50,'Exc.','Rotation',0,'FontSize',12,'FontName','Arial')
aesthetics

subplot(1,2,2)
imagesc(am(length(am):-1:1,:));
aesthetics
cb= colorbar;
aesthetics
cb.Location = 'Southoutside';
cb.Label.String = 'Correlation'
ylabel('Neuron');
xlabel('Neuron');
yts = yticks;
yticks([1 200:200:800])
yticklabels(flip(yts));
ax=gca;
ax.YLabel.Color='w';
ax.YTickLabel = [];
f=gcf;
f.Position = [700   435   596   350];
text(length(v)+20,length(v)- mean([length(v) Ne]),'Inh.','Rotation',0,'FontSize',12,'FontName','Arial')
text(length(v)+20,length(v)- mean([0 Ne]),'Exc.','Rotation',0,'FontSize',12,'FontName','Arial')
text(mean([length(v) Ne])-50,length(v)-length(v)-50,'Inh.','Rotation',0,'FontSize',12,'FontName','Arial')
text(mean([0 Ne]),length(v)-length(v)-50,'Exc.','Rotation',0,'FontSize',12,'FontName','Arial')
hold on
plot(1:length(v),Ni*ones(length(v),1),'-k','LineWidth',1.5);
plot(Ne*ones(length(v),1),1:length(v),'-k','LineWidth',1.5);
% gap = 6;
% 
% [xy1] = [length(v)-Ni, 1];
% [wh] = [Ni, Ni-gap];
% r1 = rectangle('Position',[xy1 wh],'EdgeColor','b','LineWidth',1.5);
% 
% [xy1] = [length(v)-Ni, Ni+gap];
% [wh] = [Ni, Ne-20];
% r2 = rectangle('Position',[xy1 wh],'EdgeColor','r','LineWidth',1.5);

%% 
meanIIfc = mean(mean(am(1000:-1:801,1000:-1:801)));
meanEEfc = mean(mean(am(1:800,1:800)));
    am1=am;
    am1(1000:-1:801,1000:-1:801) = NaN;
    am1(1:800,1:800) = NaN;
meanIEfc = mean(nanmean(am1));

sdIIfc = std(mean(am(1000:-1:801,1000:-1:801)));
sdEEfc = std(mean(am(1:800,1:800)));
sdIEfc = std(nanmean(am1));

round([meanIIfc sdIIfc; meanIEfc sdIEfc; meanEEfc sdEEfc],2)



