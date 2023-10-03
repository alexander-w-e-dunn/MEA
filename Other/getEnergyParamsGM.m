%  get mea params for gen model
clear all
load('MEAenergy28.mat')
load('MEAparameterspace.mat')
a = min(MEAenergy28,[],3);

for rule = 1 : size(a,2)
    for culture = 1 : size(a,1)
        sim = min (find( MEAenergy28(culture,rule,:) == a(culture,rule) ));
        param1(culture, rule) = MEAparameterspace(sim,1);
        param2(culture, rule) = MEAparameterspace(sim,2);
        clear sim
    end
end

param1mean = mean(param1);
param1mean = round(param1mean,2)';
param1sem  = std(param1) ./ (sqrt(size(a,1)));
param1sem  = round(param1sem,2)';

param2mean = mean(param2);
param2mean = round(param2mean,2)';
param2sem  = std(param2) ./ (sqrt(size(a,1)));
param2sem  = round(param2sem,2)';