%% D_o_NMC622_TK_CACold
% NMC Diffusion Coefficient Function
    % X = intercalation fraction                                        [-]
    % T = Temperature                                                   [K]
    %
    % Data taken from Todd Kingston's GITT work
function [D_o] = D_o_NMC622_TK_CACold(X , T)
    x_vec = [0.851700000000000	0.559700000000000	0.523200000000000	0.486700000000000	0.450200000000000	0.413700000000000	0.377200000000000	0.340700000000000	0.121700000000000];
    D_vec = [1.0881e-11,        1.0881e-11,         1.0687e-11,         1.4975e-11,         2.0998e-11,         2.6632e-11,         3.2084e-11,         3.85e-11,           3.85e-11         ]*1e-4;
    D_o = interp1(x_vec,D_vec,X);
end