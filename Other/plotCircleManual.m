% Manual circle plot
% AWE Dunn Cambrudge University 2020
%
% to add: 
% option to plot by module affiliation; requires inputting the module
% affiliation vector and using this to reorder the matrix then colouring
% the nodes according to module affiliation
% or colour according to manual input (e.g. cell type)
% 
% add lines linking nodes; add arrows if directional matrix
% line thickness according to connection weight
% take node labels as an input (if no labels provided don't display node
% labels

addpath(genpath('C:\Users\alexd\OneDrive - University Of Cambridge\Cam\PhD\Project\Data&CodeBackup\Scripts'))
clear all
close all
radius = 20;
centreX = 0;
centreY = 0;
numPoints = 60;
angle = 360/numPoints;
f1 = figure; scatter(centreX,centreY,'+k'); hold on
f1.Position = [680   464   560   496];
textOffset = -0.75;
for i = 1 : numPoints
   x = radius * cosd(angle*i) + centreX; %cosd gives degrees, cos function guves radians
   y = radius * sind(angle*i) + centreY;
   xtext =  (radius/10+radius) * cosd(angle*i) + centreX + textOffset;
   ytext =  (radius/10+radius) * sind(angle*i) + centreY;
   scatter(x,y,50,'filled');
   text(xtext,ytext,num2str(i),'FontSize',12);
   xy(i,[1 2]) = [x y];
end
hold off
aesthetics; box off; axis off
% figure; scatter(xy(:,1),xy(:,2)); %xlim([-1.5 1.5]); ylim([-1.5 1.5]);
%% draw curved lines between points
% 
simMat = rand(60,60);
simMat(simMat >=0.9) = 1;
simMat(simMat ~= 1) = 0;

%%
% x = [0 0 0 0 0];
% y = [-1844 -952 0 952 1844];
% z = [2000 1446 1238 1446 2000];
% 
% n = 50;
% x1 = zeros(1,n);
% y1 = linspace(min(y),max(y),n);
% z1 = interp1(y,z,y1,'spline');
% z2 = interp1(y,z,y1,'pchip');
% 
% plot3(x,y,z,'b',x1,y1,z1,'r',x1,y1,z2,'g');
% view(90,0);