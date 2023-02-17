%% Convert No Noise Step Continuous to Discrete
clear all; close all; clc;
% Desired Sampling Times
N_t_s   = 25;   % **Keep this fix for now
T_s_min = -1; % T_s = 10^(T_s_min), **Keep this fix for now
T_s_max =  1; % T_s = 10^(T_s_max), **Keep this fix for now
Ts_vec = logspace(T_s_min,T_s_max,N_t_s);

% SOC_vec = 0:1:100;

% Remove these from SOC_vec
idx_SOC = [1:80 82:90 92:94 96:97 99:101];
SOC_vec = SOC_vec(idx_SOC);

SIM = load('F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_SS_EIS_SOC50.mat','SIM');
OutputMatrix = SIM.SIM.OutputMatrix;

P = load('F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_SS_EIS_SOC50.mat','P');
P = P.P;

N_samples = 600; 
%%%% SOC_vec = 97;
for SS = 1:length(SOC_vec)
    SOC = SOC_vec(SS);
    % Create filename to load
    filename = ['F:\TylerFiles\GitHubRepos\BatteryModel\Model\Results\ObservabilityTest\ObservabilityTest_KPCont_Relax_StepSOC' num2str(SOC) '.mat'];
    data = load(filename);
    
%         % Plot Data
%         z = OutputMatrix*data.SOLN.y;
%         figure
%         plot(data.SOLN.x,z(1,:),'-o')

    for TT = 1:length(Ts_vec)
%         TT
        Ts = Ts_vec(TT);
        % Get just the CC part and then add on the rest later
        t_final_Ts = Ts*N_samples;
        t_CC_vec = 0:Ts:t_final_Ts;

        x = deval(data.SOLN,t_CC_vec);

        z = OutputMatrix * x;

        % Add Rest Data
        t_sim_vec = [0 Ts (t_CC_vec+2*Ts)];
        x = [x(:,1), x(:,1) , x];
        z = [z(:,1), z(:,1) , z];

%         % Plot Data
%         figure
%         plot(t_sim_vec,z(1,:),'-o')

        save_filepath = 'F:\TylerFiles\GitHubRepos\BatteryModel\Model\KalmanFilter\DATA\ODE_Step';
        save_filename = ['Step_SOC' num2str(SOC) '_Ts' num2str(Ts) '.mat'];
        save([save_filepath filesep save_filename],'t_sim_vec','x','z','P','OutputMatrix','SOC','Ts')

        clear x z 
    end
end




