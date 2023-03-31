clear all; close all; clc;
% I supressed KBSOC and clear all; close all; clc; in CreateProject
% Then call this to loop through each SOC for a Ts
%% Create a batch of KnownProfile
% N_SOC = 101; % This will be weird to code
% SOC_min = 0;
% SOC_max = 100;
% 
% % SOC_vec = linspace(SOC_min,SOC_max,N_SOC);
% % SOC_vec = SOC_min:5:SOC_max;
% SOC_vec = [50];
% 
% for i = 1:length(SOC_vec)
% %     i
%     % KBSOC    = SOC_vec(i);
%     PRBS_SOC = SOC_vec(i);
%     % EIS_SOC  = SOC_vec(i);
%     % SS_SOC   = SOC_vec(i);
%     CreateProject
% end

%% Loop over Sampling Rates
N_t_s = 25;
T_s_min = -1; % T_s = 10^(T_s_min)
T_s_max =  1; % T_s = 10^(T_s_max)

T_s_vec = logspace(T_s_min,T_s_max,N_t_s);
% T_s_vec = T_s_vec(2);

for i = 1:length(T_s_vec)
    % i
    PRBS_Tswitch = 10*T_s_vec(i);
    CreateProject
end