%% Convert Data from time to SOC
clear all; close all; clc;
%%
filename = 'F:\TylerFiles\GitHubRepos\BatteryModelingExtras\OLD Results\Lui_Wiley_Paper\Nondestructive_Exp_Better5C\Nondestructive_Exp_Better5C_Polar_0.20C_C.mat';
load(filename)

plot(sorted_data(:,1),sorted_data(:,2))

%% SOC Calculation
% % Cell Capacity
% Cap = -SIM.A_c/3600 * cumtrapz( t_soln , i_user );   
% 
% % SOC
% SOC = Cap/SIM.Cell_Cap + SIM.SOC_start;
t_max = 3600*5;
SOC = sorted_data(:,1)/t_max * 100;
cell_voltage = sorted_data(:,2);

V_lim = [3.4 , 4.2];

xdata = linspace(0,100,101);
ydata = zeros(size(xdata));
ydata(2:end-1) = interp1(SOC,cell_voltage,xdata(2:end-1));

ydata(1)   = V_lim(1);
ydata(end) = V_lim(2);

%% Corrent NaNs
idx = find(isnan(ydata));
ydata(idx) = interp1([xdata(idx(1)-1) , xdata(end)],[ydata(idx(1)-1) , ydata(end)],xdata(idx));


%% Plot Check
plot(xdata,ydata);

%% Save Data
save_filename = 'WileyData.mat';
save(save_filename,'xdata','ydata')