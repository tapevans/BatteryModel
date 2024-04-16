%% D_o_NMC622_TK_ANCold
% NMC Diffusion Coefficient Function
    % X = intercalation fraction                                        [-]
    % T = Temperature                                                   [K]
    %
    % Data taken from Todd Kingston's GITT work
function [D_o] = D_o_NMC622_TK_ANCold(X , T)
    x_vec = [0.851700000000000	0.559700000000000	0.523200000000000	0.486700000000000	0.450200000000000	0.413700000000000	0.377200000000000	0.340700000000000	0.121700000000000];
    D_vec = [1.1579e-11,        1.1579e-11,         1.4174e-11,         1.9891e-11,         2.6474e-11,         3.3038e-11,         3.9507e-11,         4.21e-11,           4.21e-11         ]*1e-4;
    D_o = interp1(x_vec,D_vec,X);
end