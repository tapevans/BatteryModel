%% Test Surf
close all;
% vec = [0:1:5    ;
%       (0:1:5)*2;
%       (0:1:5)*3];
% 
% Z = myHankel(vec);

Q_vec = 1:9;
R_vec = [10,20,30];

[X , Y] = meshgrid(R_vec , Q_vec);
Z = X.*Y;
surf(X,Y,Z)
xlabel('R(x)')
ylabel('Q(y)')