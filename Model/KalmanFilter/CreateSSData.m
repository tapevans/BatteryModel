%% Generate SS systems and save to their own file

clear all; close all; clc

%%
% Desired Sampling Times
N_t_s   = 25;   % **Keep this fix for now
T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
Ts_vec = logspace(T_s_min,T_s_max,N_t_s);

% Split up like this for Github pushing
SOC_vec = 0:1:25; 
% SOC_vec = 26:1:49;
% SOC_vec = 51:1:75;
% SOC_vec = 76:1:100;

% SOC_vec = 50;

for SS = 1:length(SOC_vec)
    SOC = SOC_vec(SS);
    EIS_data = load(['F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_SS_EIS_SOC' num2str(SOC) '.mat']);
    % Create CT SS
        sys = dss(EIS_data.A,EIS_data.B,EIS_data.C,EIS_data.D,EIS_data.SIM.M);
        sys_CT = ss(sys,'explicit');
        IC = EIS_data.SIM.SV_IC;

        save_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\SS_CT';
        save_filename = ['SS_CT_SOC' num2str(SOC) '.mat'];
        save([save_filepath filesep save_filename],'sys_CT','IC')

    for TT = 1:length(Ts_vec)
        Ts = Ts_vec(TT);
        % Create DT SS for each Ts
            sys_DT = c2d(sys_CT,Ts);

        save_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\SS_DT';
        save_filename = ['SS_DT_SOC' num2str(SOC) '_Ts' num2str(Ts) '.mat'];
        save([save_filepath filesep save_filename],'sys_DT','Ts','IC')
    end
end