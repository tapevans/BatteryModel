%% TDF_EC_DMC_1_1_LiPF6
% Inputs:
%   * Ce ~ Electrolyte salt concentration [kmol / m^3]
%   * T  ~ Control volume temperature [K]
%
% Outputs:
%   * TDF ~ Thermodynamic Factor [-]
%
% DOI 10.1149/2.0571912jes

function TDF = activity_Landesfeind(Ce,T)
    p_1 = -5.58E+00;
    p_2 =  7.17E+00;
    p_3 =  3.80E-02;
    p_4 =  1.91E+00;
    p_5 = -6.65E-02;
    p_6 = -5.08E-05;
    p_7 =  1.10E-01;
    p_8 = -6.10E-03;
    p_9 =  1.51E-04;
    TDF = p_1 ...
        + p_2 * Ce ...
        + p_3 * T ...
        + p_4 * Ce.^2 ...
        + p_5 * Ce .* T ...
        + p_6 * T.^2 ...
        + p_7 * Ce.^3 ...
        + p_8 * Ce.^2 .* T ...
        + p_9 * Ce .* T.^2;
end