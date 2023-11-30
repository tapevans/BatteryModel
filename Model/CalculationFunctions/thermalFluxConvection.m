%% thermalFluxConvection
%
%

function [q_conv] = thermalFluxConvection( T, SIM__h, SIM__T_inf)
    q_conv = SIM__h * (T - SIM__T_inf);
    % Convection is written positive out of the system.
    % Therefore heat flux is positive when the system temperature is
    % greater than the ambient temperature.
end