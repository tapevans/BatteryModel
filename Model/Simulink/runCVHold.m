function [time_CV, voltage_CV, current_CV] = runCVHold(AN,CA,SEP,EL,SIM,CONS,P,N,FLAG,PROPS,SV_IC,i_user_IC)
%% Use this function to run a CV hold section of a simulation

%% Change Working Directory
oldFolder = cd(pwd);


%% Parse Incoming Data
% Run Time

% Initial Conditions
% SV
% i_user


%% Constants
% t_r



%% Create Plant for Tuning PID




%% Tune PID



%% Calculate Parameters for Battery System



%% Function Handle Conversion



%% Run Simulation




%% Remove Data below C/20




%% Post-Processing



%% Return to Old Working Directory
cd(oldFolder);