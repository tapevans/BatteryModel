%% D_el_EC_DMC_1_1_LiPF6
% Inputs:
%   * Ce ~ Electrolyte salt concentration [kmol / m^3]
%   * T  ~ Control volume temperature [K]
%
% Outputs:
%   * D_el ~ binary diffusion coefficient [m^2 / s]
%
% DOI 10.1149/2.0571912jes

function D_el = D_oLiion_Landesfeind_half(Ce,T)
    p_1 =  1.47E+03;
    p_2 =  1.33E+00;
    p_3 = -1.69E+03;
    p_4 = -5.63E+02;
    D_el = (p_1 * exp(p_2 * Ce) .* exp(p_3./T) .* exp(p_4 * Ce ./ T) * 1e-6)/1000;
    % D_el = D_el * 1e6 * 1000; %makes output match the graph in paper
    D_el = D_el * 0.5;
end