%% Symbolic Hankel
clear all; close all; clc;

%% Parameters
% syms A B C D x_0 y_0 u_0
% x = sym('x_%d',[3,1]);
% u = sym('u_%d',[3,1]);
% y = sym('y_%d',[3,1]);
% 
% x = [x_0 ; x]
% u = [u_0 ; u]
% y = [y_0 ; y]

syms x_0 x_1 x_2 x_3 [2 1] matrix
syms y_0 y_1 y_2 y_3 [3 1] matrix
syms u_0 u_1 u_2 u_3 [1 1] matrix
syms A [2 2] matrix
syms B [2 1] matrix
syms C [3 2] matrix
syms D [3 1] matrix

x = [x_0, x_1, x_2, x_3];
y = [y_0, y_1, y_2, y_3];
u = [u_0, u_1, u_2, u_3];

%% Define Variables
%for i = 1:length(x)-1
% for i = 1:3
%     x(:,i+1) = A*x(:,i) + B*u(:,i) ;
% end
x_1 = A*x_0 + B*u_0;
x_2 = A*x_1 + B*u_1;
x_3 = A*x_2 + B*u_2;

% for i = 1:length(y)
% for i = 1:4
%     y(:,i) = C*x(:,i) + D*u(:,i) ;
% end

y_0 = C*x_0 + D*u_0;
y_1 = C*x_1 + D*u_1;
y_2 = C*x_2 + D*u_2;
y_3 = C*x_3 + D*u_3;

y_0
y_1
y_2
y_3

%% Insert Assumptions
% x_0 = 0;
% u = 1 when k  = 0
% u = 0 when k != 0

% u_num = zeros(size(u)); u_num(1) = 1

% x_new_new = subs(x_new,[u_0,u_1,u_2,u_3],[1,0,0,0])

% y
% y_new = subs(y,u,u_num)
% y_new_new = subs(y_new,x_0,0)

x_0_new = subs(x_0,x_0,[0;0]);
x_1_new = subs(x_1,x_0,[0;0]);
x_2_new = subs(x_2,x_0,[0;0]);
x_3_new = subs(x_3,x_0,[0;0]);
x_1_new_new = subs(x_1_new,{u_0,u_1,u_2,u_3},{1,0,0,0});
x_2_new_new = subs(x_2_new,{u_0,u_1,u_2,u_3},{1,0,0,0});
x_3_new_new = subs(x_3_new,{u_0,u_1,u_2,u_3},{1,0,0,0});

y_0_new = subs(y_0,x_0,[0;0]);
y_1_new = subs(y_1,x_0,[0;0]);
y_2_new = subs(y_2,x_0,[0;0]);
y_3_new = subs(y_3,x_0,[0;0]);
y_0_new_new = subs(y_0_new,{u_0,u_1,u_2,u_3},{1,0,0,0});
y_1_new_new = subs(y_1_new,{u_0,u_1,u_2,u_3},{1,0,0,0});
y_2_new_new = subs(y_2_new,{u_0,u_1,u_2,u_3},{1,0,0,0});
y_3_new_new = subs(y_3_new,{u_0,u_1,u_2,u_3},{1,0,0,0});


%% Define Hankel
O = [C;C*A]
C = [B,A*B]

H = O*C

H_Mark = [y_1_new_new, y_2_new_new
          y_2_new_new, y_3_new_new]

%% Test
% A = sym('A_%d%d',[2,2]);
% B = sym('B_%d%d',[2,1]);
% C = sym('C_%d%d',[3,2]);
% D = sym('D_%d%d',[3,1]);
% 
% A
% B
% C
% D
% 
% C*B
% try 
%     B*C
% catch ME
%     %rethrow(ME)
%     disp('Error Occurred')
% end
% 
% syms A [2 2] matrix
% syms B [2 1] matrix
% syms C [3 2] matrix
% syms D [3 1] matrix
% A
% B
% C
% D
% 
% C*B
% try 
%     B*C
% catch ME
%     rethrow(ME)
%     %disp('Error Occurred')
% end