%% thermalFluxConvection
%
%

function [q_conv] = thermalFluxConvection( SV , SIM , P )
    T      = SV( P.T , :);
    q_conv = SIM.h * (T - SIM.T_inf);
    % Convection is written positive out of the system.
    % Therefore heat flux is positive when the system temperature is
    % greater than the ambient temperature.
end