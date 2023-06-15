%% Test Inverse CPCT
clear all; close all; clc;

%%

N_outputs = 6;
N_states  = 6; 
N_samples = 10000;

% C = 1: N_outputs*N_states; % Using this results in a rank deficient matrix
C = randn(N_outputs*N_states,1);
C = reshape(C,N_outputs,[]);

sigma = logspace(-2,-1,N_states);

P = diag(sigma.^2);

D = C';
A = C'*getPinv(C*C');
B = getPinv(D'*D)*D';

% C\C
% C*A
% B*C'
% 
% A'*C'
% D'*B'
% 
% % getPinv(C)
% % (getPinv(C)*C)'*getPinv(C)*C
% 
% % C*getPinv(C)
% % getPinv(C')*C'
% 
% C*P*C' - C*P*D
% D'*P'*C' - C*P*D

% inv(C)*C
% C'*inv(C')

H = C*P*C';
P_rev = inv(C)*H*inv(C');
diff_P = P - P_rev


% % C*A*P*B*D 
% % 
% % C*A*P*B*D - P
% % 
% % D'*B'*P'*A'*C'
% % 
% % D'*B'*P'*A'*C' - P'
% % 
% % BDT = (B*D)';
% % CAT = (C*A)';








%% Helper Function
function [X_inv] = getPinv(X)
    [U,S,V] = svd(X);
    S_inv = pinv(S);
    X_inv = V * S_inv * U';
end