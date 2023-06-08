%% Testing Wrong IC Part 2
clear all; close all; clc

%% Subdirecties to Include
    % Script's filepath
    [current_file_path,~,~] = fileparts(mfilename('fullpath'));

    % Ensure this is the current directory, important for pwd
    cd(current_file_path)
    
    % Include all folders
    addpath(genpath(current_file_path));
    % genpath creates a string with all folders and subfolders in a given directory
    % addpath then adds all of them to the current workspace


%% Parameters
    filename = 'ROM_SS_SOC50_Ts1.mat';

    
%% Load SS
    load(filename)


%% System Variables
    sys = sys_HK{end}; % This is hardcoded for combined system for right now

    A_DT    = sys.A;
    B_DT    = sys.B;
    C_DT    = sys.C;
    C_DT_CV = sys.C(1,:);

    [N.measur, N.states] = size(C_DT_CV);
    [N.measur_all, ~]    = size(C_DT);


%% Desired Deviation
    sigma = 3e-1;
    y_des = sigma * ones(N.measur_all,1);
    
    % Try using inverse
    x_hat = C_DT \ y_des;

    % Calculate P_k|k
    error = zeros(size(x_hat)) - x_hat;

    P_kk = nan(length(x_hat));
    for i = 1:length(x_hat)
        for j = 1:length(x_hat)
            [P_kk(i,j)] = calcPopCovar(error(i,:), error(j,:));
        end
    end

    % Calculate Output Covar
    CPCT = C_DT * P_kk * C_DT';
    sigma_calc = sqrt(CPCT(1,1));




    
    


%% 
function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end



