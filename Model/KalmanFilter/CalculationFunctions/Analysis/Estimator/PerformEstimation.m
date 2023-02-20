%% Estimator

function [ESTIMATOR,RESULTS] = PerformEstimation(plant_filename,SIM,FLAG,N,P,RESULTS)
%% Load Plant Data
Plant = load(plant_filename);

% Plant Data
u = Plant.Plant_Data.Fs.Fs_soln;
z = Plant.Plant_Data.zN.z_soln;


%% Get Estimator Model
%  1) Matlab SS_DT
%  2) Ho-Kalman

% Overall System
    old_C_mode = FLAG.C_mode;
    FLAG.C_mode = 5;
    switch FLAG.EstimatorModel
        case 1
            [~ , sys_DT] = getSS_System(SIM,N,P,FLAG);
            est_sys = sys_DT;
        case 2
            [est_sys] = getHoKalmanROM(SIM,N,P,FLAG);
    end
    FLAG.C_mode = old_C_mode;
    C_des_ROM = est_sys.C;
    ESTIMATOR.C_des_ROM = C_des_ROM;
    ESTIMATOR.sys = est_sys;

% Reduce the outputs to C_mode
    switch FLAG.C_mode
        case N.states
            C_m = est_sys.C;
            D_m = est_sys.D;
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
        otherwise
            C_m = est_sys.C(FLAG.C_mode,:);
            D_m = est_sys.D(FLAG.C_mode,:);
            est_sys = ss(est_sys.A,est_sys.B,C_m,D_m,SIM.Ts);
    end

% Initial Conditions
    y_0 = SIM.x_0;
    x_hat_0 = C_des_ROM\y_0;
%     x_hat_0(3) = 5.25;
%     x_hat_0(4) = -0.25;


%% Do Asymptotic Calcs
[ESTIMATOR.K_infty, ESTIMATOR.P_infty] = AsymptoticPreCalcs(FLAG,SIM,est_sys);


%% Do Asymptotic Estimation
[x_hat_asy] = AsyEstimator(est_sys,FLAG,SIM,x_hat_0,u,z);
y_hat_asy = est_sys.C * x_hat_asy;
y_hat_asy_ALL = C_des_ROM * x_hat_asy;

RESULTS.EST.ASY.t_soln = Plant.Plant_Data.Fs.t_soln;
RESULTS.EST.ASY.x_soln = x_hat_asy;
RESULTS.EST.ASY.z_soln = y_hat_asy;
RESULTS.EST.ASY.z_soln_ALL = y_hat_asy_ALL;


%% Do Variable Estimation
[x_hat_var, ESTIMATOR.K_k, ESTIMATOR.P_k_pre] = VarEstimator(est_sys,FLAG,SIM,x_hat_0,u,z);
y_hat_var = est_sys.C * x_hat_var;
y_hat_var_ALL = C_des_ROM * x_hat_var;

RESULTS.EST.VAR.t_soln = Plant.Plant_Data.Fs.t_soln;
RESULTS.EST.VAR.x_soln = x_hat_var;
RESULTS.EST.VAR.z_soln = y_hat_var;
RESULTS.EST.VAR.z_soln_ALL = y_hat_var_ALL;


%% Convert Plant true states into ROM states
[U, S, V]   = svd(C_des_ROM);
threshold   = 1e-7;
S_cross     = pinv(S,threshold);
C_des_cross = V*S_cross*U';

x_plant_ROM = C_des_cross * Plant.Plant_Data.zNA.z_soln;

RESULTS.EST.PLANT.t_soln = Plant.Plant_Data.Fs.t_soln;
RESULTS.EST.PLANT.x_soln = x_plant_ROM;
RESULTS.EST.PLANT.z_soln = Plant.Plant_Data.zN.z_soln;


end