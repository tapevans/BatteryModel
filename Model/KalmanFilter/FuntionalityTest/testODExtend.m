%% Testing odextend
% Results:
% For this specific example, it didn't work. I am going to test a different
% ODE


clear all; close all; clc;

%%
% https://www.mathworks.com/help/matlab/ref/odextend.html
tspan = [0 20];
y0    = [2 0];
sol   = ode45(@vdp1,tspan,y0);
% sol = struct with fields:
%      solver: 'ode45'
%     extdata: [1x1 struct]
%           x: [0 1.0048e-04 6.0285e-04 0.0031 0.0157 0.0785 0.2844 0.5407 0.8788 1.4032 1.8905 2.3778 2.7795 3.1285 3.4093 3.6657 3.9275 4.2944 4.9013 5.3506 5.7998 6.2075 6.5387 6.7519 6.9652 7.2247 7.5719 8.1226 8.6122 9.1017 9.5054 9.8402 10.1157 ... ]
%           y: [2x60 double]
%       stats: [1x1 struct]
%       idata: [1x1 struct]

% Use linspace to generate 250 points in the interval [0 20]. Evaluate the solution at these points using deval.
x = linspace(0,20,250);
y = deval(sol,x);

% Plot the first component of the solution.
plot(x , y(1,:))



% Extend the solution to tf=35 using odextend and add the result to the original plot.
sol_ext = odextend(sol,@vdp1,35);
x_ext   = linspace(20,35,350);
y_ext   = deval(sol_ext,x_ext);
hold on
plot(x_ext , y_ext(1,:) , 'r')


%% Modification
% I don't want sol_new to store all of the previous data as well

% Modify sol
n_points = length(sol.x);
% y_IDX    = [1:n_points];
% y_IDX    = [1, n_points-10:n_points];
y_IDX    = [n_points-1:n_points];

sol_mod.solver  = sol.solver;
sol_mod.extdata = sol.extdata;
sol_mod.x       = sol.x(1,y_IDX);
sol_mod.y       = sol.y(:,y_IDX);
sol_mod.stats   = sol.stats;
sol_mod.idata   = sol.idata;

y_IC = sol.y(:,end);

% Extend the solution to tf=35 using odextend and add the result to the original plot.
sol_mod_ext = odextend(sol_mod,@vdp1,35,y_IC);
x_mod_ext   = linspace(20,35,350);
y_mod_ext   = deval(sol_mod_ext,x_mod_ext);
hold on
plot(x_mod_ext , y_mod_ext(1,:) , 'ko')







