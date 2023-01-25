%% Testing Simulink Noise generation
%
%
clc;
% close all;
%% 
N_inputs = 1;
Q = 1e-4;
seed = 1;
rng(seed);
w_k = (chol(Q,'lower')*randn(N_inputs,201))';
% w_k = (chol(Q,'lower')*randn(N_inputs,1e6))';

figure
hold on
plot(w_k,'k-')
plot(out.w_k,'bo') % have to run the simulink model to get this data
xlim([0,200])


% It seems that simulink and Matlab don't produce the same random number
% vector