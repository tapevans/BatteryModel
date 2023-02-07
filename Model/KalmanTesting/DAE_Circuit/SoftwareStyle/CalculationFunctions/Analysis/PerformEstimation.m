%% Estimator

function [x_hat_asy, x_hat_var, y_hat_asy, y_hat_var, P_infty, K_infty, K_k, P_k_pre] = PerformEstimation(filename)
%% Load Plant Data
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\Plant3_IQ_Step_Q1_R1_Ts1e-06_CAll.mat';
Plant = load(filename);

% Plant Data
u = Plant.u;
z = Plant.z;
sys = Plant.sys_DT3;
y_0 = Plant.SIM.x_0_3;
x_hat_0 = Plant.sys_CT3.C\y_0;
% x_hat_0 = [0;0;0];

% sys = Plant.sys_DT5;
FLAG = Plant.FLAG;
SIM  = Plant.SIM;

%% Do Asymptotic Calcs
[K_infty, P_infty] = AsymptoticPreCalcs(FLAG,SIM,sys);

%% Do Asymptotic Estimation
[x_hat_asy] = AsyEstimator(sys,FLAG,SIM,x_hat_0,u,z);
y_hat_asy = sys.C * x_hat_asy;

%% Do Variable Estimation
[x_hat_var, K_k, P_k_pre] = VarEstimator(sys,FLAG,SIM,x_hat_0,u,z);
y_hat_var = sys.C * x_hat_var;

end