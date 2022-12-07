clear all; close all; clc;
% I supressed KBSOC and clear all; close all; clc; in CreateProject
% Then call this to loop through each SOC for a Ts
%% Create a batch of KnownProfile
% N_t_s = 25;
% T_s_min = -1; % T_s = 10^(T_s_min)
% T_s_max =  1; % T_s = 10^(T_s_max)

N_SOC = 101; % This will be weird to code
SOC_min = 0;
SOC_max = 100;

% T_s_vec = logspace(T_s_min,T_s_max,N_t_s);
SOC_vec = linspace(SOC_min,SOC_max,N_SOC);

for i = 1:length(SOC_vec)
%     i
    KBSOC = SOC_vec(i);
    CreateProject
end