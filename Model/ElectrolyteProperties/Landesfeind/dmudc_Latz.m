%% dmudc
% Inputs:
%   * Ce  ~ Electrolyte salt concentration [kmol / m^3]
%   * T   ~ Control volume temperature [K]
%   * tf  ~ Transference Number [-]
%   * TDF ~ Thermodynamic Factor [-]
%
% Outputs:
%   * dmudc ~ Latz thermodynamic factor [J m^3 / kmol^2]

function dmudc = dmudc_Latz( Ce , T, tf , TDF)
    R     = 8314.472; % [J kmol^-1 K^-1]
    dmudc = 2*R.*T .* TDF .* (1 - tf) ./ (tf .* Ce);
end