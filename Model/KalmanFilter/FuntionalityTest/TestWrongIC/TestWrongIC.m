%% Test Setting the initial wrong conditions
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


%% Test Calcs
    % var = 1e-4;
    % NN = 1e6;
    % noise        = sqrt(var)*randn(1,NN);
    % avg          = mean(noise)
    % diff         = noise - avg;
    % squ_diff     = diff.^2;
    % sum_squ_diff = sum(squ_diff);
    % var_calc     = sum_squ_diff/NN
    % std          = sqrt(var)
    % covar        = calcPopCovar(noise, noise)

    y_expected_single = (1:N.measur_all)'; % Expected Value
    y_expected_single = [3.749 ; 0.101 ; 0 ; 1 ; 14.80 ; 0];

    var = 1e-4;
    % y_estimated_single = y_expected_single*sqrt(offset) + y_expected_single;
    y_estimated_single = sqrt(var) + y_expected_single;
    numberSamples = 1;
    y_estimated = [];
    y_expected  = [];
    for i = 1:numberSamples
        y_estimated = [y_estimated , y_estimated_single];
        y_expected  = [y_expected  , y_expected_single];
    end
    
    % Calculate Output Covariance
    error_outputs = y_expected - y_estimated;
    
    covar_outputs = nan(N.measur_all);
    for i = 1:N.measur_all
        for j = 1:N.measur_all
            [covar_outputs(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
        end
    end

    % Calculate State Covariance
    [U,S,V] = svd(C_DT);
    S_inv = pinv(S);
    C_DT_inv = V*S_inv*U';

    [U,S,V] = svd(C_DT');
    S_inv = pinv(S);
    C_DT_tpose_inv =  V*S_inv*U';

    P_kk = C_DT_inv * covar_outputs * C_DT_tpose_inv;
    
    % State Variance
        var_state = diag(P_kk);

    % State STD
        STD_state = (var_state).^0.5;

    % State Initial Deviation
        StateIC_offset = 3*STD_state;

    % Initial States
        state_IC     =    StateIC_offset;    
        % state_IC     =    StateIC_offset.^(1/10);
        % state_IC(2)  = -1*state_IC(2);
        % state_IC(4)  = -1*state_IC(4);
        % state_IC(6)  = -1*state_IC(6);
        % state_IC(8)  = -1*state_IC(8);
        % state_IC(10) = -1*state_IC(10);


%% Test Initial Conditions Output Covariance
    % Wrong Outputs
    y_IC_wrong = C_DT * state_IC + y_expected(:,1);

    % Calculate Output Covariance
    error_outputs_test = y_expected - y_IC_wrong;
    
    covar_outputs_test = nan(N.measur_all);
    for i = 1:N.measur_all
        for j = 1:N.measur_all
            [covar_outputs_test(i,j)] = calcPopCovar(error_outputs_test(i,:), error_outputs_test(j,:));
        end
    end

    howWrong = covar_outputs_test./covar_outputs;

%% Reverse the analysis
% Set wrong states
    x_expected_single = zeros(10,1);
    var = 1e-4;
    x_estimated_single = sqrt(var) + x_expected_single;
    numberSamples = 5;
    x_estimated = [];
    x_expected  = [];
    for i = 1:numberSamples
        x_estimated = [x_estimated , x_estimated_single];
        x_expected  = [x_expected  , x_expected_single];
    end
    
% Calculate State Error Covariance
    error_outputs = x_expected - x_estimated;
    
    P_kk = nan(length(x_expected_single));
    for i = 1:length(x_expected_single)
        for j = 1:length(x_expected_single)
            [P_kk(i,j)] = calcPopCovar(error_outputs(i,:), error_outputs(j,:));
        end
    end

% Check Calculation
    % y_tilde = ( C_DT*(x_expected - x_estimated) ) ;
    y_expected = C_DT*x_expected;
    y_estimated = C_DT*x_estimated;

    error_outputs_test = y_expected - y_estimated;

    CPC_test = nan(N.measur_all);
    for i = 1:N.measur_all
        for j = 1:N.measur_all
            [CPC_test(i,j)] = calcPopCovar(error_outputs_test(i,:), error_outputs_test(j,:));
        end
    end

    out_error_covar = C_DT * P_kk * C_DT';

    howWrong_Part2 = CPC_test./out_error_covar;



    


%% 
function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end


