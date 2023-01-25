%% Noise Testing
clear all;
close all;
clc;

%% 
% Q = 1e-4*(eye(2));
Q = [4e-1 0
     0 1e3];
% N_length = 200;
N_length = 1e6;

%% sqrt(q)*randn
w_sqrt = sqrt(Q)*randn(length(Q),N_length);
[mu_x] = calcMean(w_sqrt(1,:));
[mu_y] = calcMean(w_sqrt(2,:));
[error_x] = calcError(w_sqrt(1,:),mu_x);
[error_y] = calcError(w_sqrt(2,:),mu_y);
[Covar_xy] = calcPopCovar(error_x, error_y);
[Covar_yx] = calcPopCovar(error_y, error_x);
[Covar_xx] = calcPopCovar(error_x, error_x);
[Covar_yy] = calcPopCovar(error_y, error_y);
Covar = [Covar_xx Covar_xy
         Covar_yx Covar_yy]

%% chol(Q)
w_k_chol = (chol(Q,'lower')*randn(length(Q),N_length));
[mu_x] = calcMean(w_k_chol(1,:));
[mu_y] = calcMean(w_k_chol(2,:));
[error_x] = calcError(w_k_chol(1,:),mu_x);
[error_y] = calcError(w_k_chol(2,:),mu_y);
[Covar_xy] = calcPopCovar(error_x, error_y);
[Covar_yx] = calcPopCovar(error_y, error_x);
[Covar_xx] = calcPopCovar(error_x, error_x);
[Covar_yy] = calcPopCovar(error_y, error_y);
Covar = [Covar_xx Covar_xy
         Covar_yx Covar_yy]


%% Run Simulink
    mdl = 'TestSimulinkNoise';
    load_system(mdl)
    in = Simulink.SimulationInput(mdl);
    in = in.setModelParameter('StartTime','0','StopTime',num2str(N_length));
    mdlWks = get_param(in,'ModelWorkspace');
%     assignin(mdlWks,'Q'      ,Q)
    assignin(mdlWks,'Q'      ,diag(Q))
    out = sim(in);

    
%% Simulink randn
w_k_rand(1,:) = out.w_k_rand(1,1,:);
% w_k_rand(2,:) = out.w_k_rand(2,2,:);
w_k_rand(2,:) = out.w_k_rand(2,1,:);
[mu_x] = calcMean(w_k_rand(1,:));
[mu_y] = calcMean(w_k_rand(2,:));
[error_x] = calcError(w_k_rand(1,:),mu_x);
[error_y] = calcError(w_k_rand(2,:),mu_y);
[Covar_xy] = calcPopCovar(error_x, error_y);
[Covar_yx] = calcPopCovar(error_y, error_x);
[Covar_xx] = calcPopCovar(error_x, error_x);
[Covar_yy] = calcPopCovar(error_y, error_y);
Covar = [Covar_xx Covar_xy
         Covar_yx Covar_yy]


%% wng

%% Simulink wng


%% Functions
function [mu] = calcMean(x)
    mu = mean(x);
end
function [error] = calcError(x,mu)
% difference between x and the expected value (mu)
    error = x-mu;
end
function [Covar] = calcPopCovar(error_x, error_y)
    Covar = (error_x*error_y')/length(error_x);
end
function [Covar] = calcSampleCovar(error_x, error_y)
    Covar = (error_x*error_y')/(length(error_x)-1);
end


%% OOOOOOOOOOOOOOOOLD
%% Test vector from wiki
% X = [5 5 5 6 6 6 6 7 7 7 5 5 5 6 6 6 6 7 7 7];
% Y = [9 9 9 8 8 8 8 8 9 9 9 9 9 8 8 8 8 8 9 9];
% [mu_x] = calcMean(X);
% [mu_y] = calcMean(Y);
% [error_x] = calcError(X,mu_x);
% [error_y] = calcError(Y,mu_y);
% [Covar_xy] = calcCovar(error_x, error_y);
% [Covar_yx] = calcCovar(error_y, error_x);
% [Covar_xx] = calcCovar(error_x, error_x);
% [Covar_yy] = calcCovar(error_y, error_y);
% 
% Covar = [Covar_xx Covar_xy
%          Covar_yx Covar_yy]
% Covar_sqrt = Covar.^(1/2)

%% Another Test
% X = [5 12 18 23 45];
% Y = [2 8  18 20 28];
% [mu_x] = calcMean(X);
% [mu_y] = calcMean(Y);
% [error_x] = calcError(X,mu_x);
% [error_y] = calcError(Y,mu_y);
% [Covar_xy] = calcPopCovar(error_x, error_y);
% [Covar_yx] = calcPopCovar(error_y, error_x);
% [Covar_xx] = calcPopCovar(error_x, error_x);
% [Covar_yy] = calcPopCovar(error_y, error_y);
% Covar = [Covar_xx Covar_xy
%          Covar_yx Covar_yy]
% 
% [Covar_xy] = calcSampleCovar(error_x, error_y);
% [Covar_yx] = calcSampleCovar(error_y, error_x);
% [Covar_xx] = calcSampleCovar(error_x, error_x);
% [Covar_yy] = calcSampleCovar(error_y, error_y);
% Covar = [Covar_xx Covar_xy
%          Covar_yx Covar_yy]
