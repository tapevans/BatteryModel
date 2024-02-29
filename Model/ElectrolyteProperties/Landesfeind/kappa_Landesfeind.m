%% kappa_EC_DMC_1_1_LiPF6
% Inputs:
%   * Ce ~ Electrolyte salt concentration [kmol / m^3]
%   * T  ~ Control volume temperature [K]
%
% Outputs:
%   * kappa ~ ionic conductivity [S/m]
%
% DOI 10.1149/2.0571912jes
function kappa = kappa_Landesfeind(Ce,T)
    p_1 =  7.98e-1;
    p_2 =  2.28e+2;
    p_3 = -1.22e+0;
    p_4 =  5.09e-1;
    p_5 = -4.00e-3;
    p_6 =  3.79e-3;
    kappa = (p_1 * (1 + (T - p_2)) .* Ce .* ( 1  +  p_3 * sqrt(Ce)  +  p_4 * ( 1 + p_5 * exp(1000./T)) .* Ce ) ./ ( 1 + Ce.^4 .* (p_6 * exp(1000./T)) ))/10;
end