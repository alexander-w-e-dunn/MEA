% 
% 
% simulate gaussian noise, how many spikes would you detect
% could try red noise as it's more biological
% 

r = normrnd(mean(filtdata),std(filtdata),length(filtdata),1);
thr = mean(r) - (3*std(r));
thr
length(find(r<thr))
figure;histogram(r)
hold on, h = histogram(filddata,'b')
hold on, h = histogram(filtdata,'b')
hold on, h = histogram(filtdata)
figure;h1 = histogram(r)
hold on, h2 = histogram(filtdata)
figure;h1 = histogram(r)
hold on, h2 = histogram(filtdata)
hold on; h2 = histogram(filtdata)
h2
h2.FaceColor = 'r'
hold on; h2 = histogram(filtdata,'FaceColor','r')
close all
hold on; h2 = histogram(filtdata,'FaceColor','r')
hold on; h2 = histogram(filtdata,'EdgeColor','r')
hold on; h2 = histogram(filtdata,'EdgeColor','r','EdgeAlpha',0.2)
figure;h1 = histogram(r)
hold on; h2 = histogram(filtdata,'EdgeColor','r','EdgeAlpha',0.2)
length(find(r<-13))
length(find(filtdata<-13))
440 / (4500000)/25000))
440 / (4500000/25000)
13/std(filtdata)
thr = mean(r) - (3*std(r));
length(find(r<thr))
length(find(filtdata<thr))
iqr(filtdata)
std(filtdata)
thr = median(r) - (3*iqr(r))
length(find(r<thr))
length(find(filtdata<thr))
set(gca,'YScale','log')
thr = median(filtdata) - (3*iqr(filtdata))
length(find(r<thr))
length(find(filtdata<thr))
thr = mean(filtdata) - (3*std(filtdata))
length(find(r<thr))
length(find(filtdata<thr))
thr = mean(filtdata) - (4*std(filtdata))
length(find(r<thr))
length(find(filtdata<thr))
mean(filtdata)
median(filtdata)