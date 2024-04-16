%% D_o_NMC622_TK
% NMC Diffusion Coefficient Function
    % X = intercalation fraction                                        [-]
    % T = Temperature                                                   [K]
    %
    % Data taken from Todd Kingston's GITT work
function [D_o] = D_o_NMC622_TK_ISO20(X , T)
    x_vec = [0.851700000000000	0.559700000000000	0.523200000000000	0.486700000000000	0.450200000000000	0.413700000000000	0.377200000000000	0.340700000000000	0.121700000000000];
    D_vec = [1.344e-11,         1.344e-11,          1.6275e-11,         2.2118e-11,         2.8958e-11,         3.572e-11,          4.2629e-11,         4.53e-11,           4.53e-11         ]*1e-4;
    D_o = interp1(x_vec,D_vec,X);
end