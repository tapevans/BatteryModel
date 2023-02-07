%% Governing Equations
% This is called by Simulink
function [dSVdt] = DAECircGovnEqn(t,SV,SIM,N,P,FLAG,V_in)
% Input Current (Step)
if FLAG.InputSignalMode == 1
    [V_in] = getInput(t,SIM);
else
    % Use the V_in from the inputs
end


% Rephrase Parameters
    C_1 = SIM.C_1;
    C_4 = SIM.C_4;
    R_1 = SIM.R_1;
    R_2 = SIM.R_2;
    R_3 = SIM.R_3;
    R_4 = SIM.R_4;

% % dSVdt Equations
%     dSVdt = zeros(N.states,1);
% 
%     dSVdt(P.V_1,1)  = - (1/R_1)*(SV(P.V_1) - SV(P.V_s)) - (1/R_2)*(SV(P.V_1) - SV(P.V_2));
%     dSVdt(P.V_2,1)  =   (1/R_2)*(SV(P.V_2) - SV(P.V_1)) + (1/R_3)*(SV(P.V_2) - SV(P.V_3));
%     dSVdt(P.V_3,1)  = - (1/R_3)*(SV(P.V_3) - SV(P.V_2)) - (1/R_4)*(SV(P.V_3) - 0        );
%     dSVdt(P.V_s,1)  = - SV(P.i_PS) - (1/R_1)*(SV(P.V_s) - SV(P.V_1));
%     dSVdt(P.i_PS,1) =   SV(P.V_s)  - V_in;

% dSVdt Equations
    dSVdt = zeros(N.states,1);

    dSVdt(P.V_s,1)     = - SV(P.i_PS) - (1/R_1)*(SV(P.deltaV));
    dSVdt(P.deltaV,1)  =   (1/R_1)*(SV(P.deltaV)) - (1/R_2)*(SV(P.V_s) - SV(P.deltaV) - SV(P.V_2));
    dSVdt(P.V_2,1)     =   (1/R_2)*(SV(P.V_2) - SV(P.V_s) + SV(P.deltaV)) + (1/R_3)*(SV(P.V_2) - SV(P.V_3));
    dSVdt(P.V_3,1)     = - (1/R_3)*(SV(P.V_3) - SV(P.V_2)) - (1/R_4)*(SV(P.V_3) - 0        );
    dSVdt(P.i_PS,1)    =   SV(P.V_s)  - V_in;

end