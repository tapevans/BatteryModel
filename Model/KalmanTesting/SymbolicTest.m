clear all; close all; clc;

syms A B C A_bar B_bar C_bar T

A = sym('A_%d%d',[3,3])
B = sym('B_%d',[3,1])
C = sym('C_%d',[1,3])

% A_bar = sym('A_bar_%d%d',[3,3])
% B_bar = sym('B_bar_%d',[3,1])
% C_bar = sym('C_bar_%d',[1,3])

T = sym('T_%d%d',[3,3])
T_inv = inv(T)

%%
A_bar = T_inv*A*T
B_bar = T_inv*B
C_bar = C*T

% eqn_C = C*T == C_bar