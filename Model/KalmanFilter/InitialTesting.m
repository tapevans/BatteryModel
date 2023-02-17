%% Testing
clear all; close all; clc;
format long
%% Load EIS File
filename_EIS = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_SS_EIS_SOC50.mat';
Data_EIS = load(filename_EIS);

%% Load Impulse Data
filename_Imp = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_KPCont_DTImpulseTs1.0SOC50.mat';
Data_Imp = load(filename_Imp);

%% Solving for how long to run a simulation to get a 1% SOC difference
% Assumes i_user is 1A/m^2
Cell_Cap = Data_EIS.SIM.Cell_Cap;
A_c      = Data_EIS.SIM.A_c;
t_1perc  = (Cell_Cap/(A_c*100))*60*60
%  t_1perc  = 0.244837153197069 hr
%  t_1perc  = 8.814137515094490e+02 s

Data_Imp.SIM.Controller_MO_File(2).C_rate;
1/Data_Imp.SIM.Controller_MO_File(2).C_rate;
(1/Data_Imp.SIM.Controller_MO_File(2).C_rate)*60*60/100

N_samples = 600;

Ts = .10;
t_step_short = Ts*N_samples
t_step_perc = t_step_short / t_1perc

Ts = 1;
t_step_mid = Ts*N_samples
t_step_perc = t_step_mid / t_1perc

Ts = 10;
t_step_long = Ts*N_samples
t_step_perc = t_step_long / t_1perc

Ts_1percent = t_1perc/N_samples

Ts = Ts_1percent;
t_step_1per = Ts*N_samples
t_step_perc = t_step_1per / t_1perc

% Take aways from this:
%  1) From the previous analysis, 600 samples is sufficient for CPC calc to
%     converge to DARE CPC.
%  2) At shorter and mid Ts, 600 samples for a step will result in less
%     than 1% change in SOC (Good if trying to stay linear)
%  3) For this setup, Ts = 1.469 results in a 1% change in SOC
%  4) At higher SOC, I should flip the simulation from charge to discharge
%     step so it is not going above 100% SOC. The largest SOC change is
%     less than 7%. I can flip everything above 90% to be discharge

%% 