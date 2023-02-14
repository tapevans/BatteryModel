%% Estimator

function [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS)
%% Load Plant Data
% filename = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanTesting\DAE_Circuit\SoftwareStyle\Results\Plant3_IQ_Step_Q1_R1_Ts1e-06_CAll.mat';
Plant = load(plant_filename);

% Plant Data
u = Plant.Plant_Data.Vs.Vs_soln;
z = Plant.Plant_Data.zN.z_soln;


%% Get Estimator Model
%  1) Matlab SS_DT3
%  2) Matlab SS_DT5
%  3) Ho-Kalman

% Overall System
    old_C_mode = FLAG.C_mode;
    FLAG.C_mode = 4;
    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT3 , ~ , ~] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT3;
        case 2
            [~ , ~ , ~ , sys_DT5] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT5;
        case 3
            % Get Ho-Kalman ROM
    end
    FLAG.C_mode = old_C_mode;
    C_des_ROM = est_sys.C;
    ESTIMATOR.C_des_ROM = C_des_ROM;
    ESTIMATOR.sys = est_sys;

% Reduce the outputs to C_mode
    switch FLAG.C_mode
        case 1
            C_m = est_sys.C(1,:);
            D_m = est_sys.D(1,:);
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
        case 2
            C_m = est_sys.C(2,:);
            D_m = est_sys.D(2,:);
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
        case 3
            C_m = est_sys.C(3,:);
            D_m = est_sys.D(3,:);
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
        case 4
            C_m = est_sys.C;
            D_m = est_sys.D;
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
    end

% Initial Conditions
    switch FLAG.SlinkModel
        case 3 % 3 state plant
            y_0 = SIM.x_0_3;
        case 5
            y_0 = SIM.x_0_5;
    end
    x_hat_0 = C_des_ROM\y_0;
    
%     x_hat_0 = [0;10];

%% Do Asymptotic Calcs
[ESTIMATOR.K_infty, ESTIMATOR.P_infty] = AsymptoticPreCalcs(FLAG,SIM,est_sys);

%% Do Asymptotic Estimation
[x_hat_asy] = AsyEstimator(est_sys,FLAG,SIM,x_hat_0,u,z);
y_hat_asy = est_sys.C * x_hat_asy;
y_hat_asy_ALL = C_des_ROM * x_hat_asy;

RESULTS.EST.ASY.t_soln = Plant.Plant_Data.Vs.t_soln;
RESULTS.EST.ASY.x_soln = x_hat_asy;
RESULTS.EST.ASY.z_soln = y_hat_asy;
RESULTS.EST.ASY.z_soln_ALL = y_hat_asy_ALL;

%% Do Variable Estimation
[x_hat_var, ESTIMATOR.K_k, ESTIMATOR.P_k_pre] = VarEstimator(est_sys,FLAG,SIM,x_hat_0,u,z);
y_hat_var = est_sys.C * x_hat_var;
y_hat_var_ALL = C_des_ROM * x_hat_var;

RESULTS.EST.VAR.t_soln = Plant.Plant_Data.Vs.t_soln;
RESULTS.EST.VAR.x_soln = x_hat_var;
RESULTS.EST.VAR.z_soln = y_hat_var;
RESULTS.EST.VAR.z_soln_ALL = y_hat_var_ALL;

%% Convert Plant 3 states into 2 (or the rank of ROM)
[U, S, V]   = svd(C_des_ROM);
threshold   = 1e-7;
S_cross     = pinv(S,threshold);
C_des_cross = V*S_cross*U';

x_plant_ROM = C_des_cross * Plant.Plant_Data.zNA.z_soln;

RESULTS.EST.PLANT.t_soln = Plant.Plant_Data.Vs.t_soln;
RESULTS.EST.PLANT.x_soln = x_plant_ROM;
RESULTS.EST.PLANT.z_soln = Plant.Plant_Data.zN.z_soln;

end