%% Transference Number
% Inputs:
%   * Ce ~ Electrolyte salt concentration [kmol / m^3]
%   * T  ~ Control volume temperature [K]
%
% Outputs:
%   * tfNum ~ Transference Number [-]
%
% DOI 10.1149/2.0571912jes
function t = transferenceNumber_Landesfeind(Ce, T)
    p_1 = -7.91E+00;
    p_2 =  2.45E-01;
    p_3 =  5.28E-02;
    p_4 =  6.98E-01;
    p_5 = -1.08E-02;
    p_6 = -8.21E-05;
    p_7 =  7.43E-04;
    p_8 = -2.22E-03;
    p_9 =  3.07E-05;
    t =   p_1 ...
        + p_2 * Ce ...
        + p_3 * T ...
        + p_4 * Ce.^2 ...
        + p_5 * Ce .* T ...
        + p_6 * T.^2 ...
        + p_7 * Ce.^3 ...
        + p_8 * Ce.^2 .* T ...
        + p_9 * Ce .* T.^2;
end