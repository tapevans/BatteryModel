%% thermalFluxConvection
%
%

function [q_conv] = thermalFluxConvection( SV , SIM , P )
    T      = SV( P.T , :);
    q_conv = SIM.h * (T - SIM.T_inf);
end