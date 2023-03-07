%% Plot Surf test instead of contour
close all;
clear all;
clc;
%%
x_vec = 1:1:5;
y_vec = 5:1:9;
[X,Y] = meshgrid(x_vec,y_vec);
Z = 1:1:25;
Z = reshape(Z,size(X));

% figure
% surf(X,Y,Z)
% 
% figure
% contourf(X,Y,Z)

figure
b = bar3(y_vec,Z,1);
set(b,'EdgeColor','none')
view(-90,90)
xlabel('x')
ylabel('y')