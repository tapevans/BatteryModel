%% Governing Equations (3 State System)
function [dSVdt] = GovnEqn_3(t,SV,SIM,N,P,FLAG,A,B,V_in)
%% Input Current (Step)
if FLAG.InputSignalMode == 1
    [V_in] = V_inCalc(t,SIM,FLAG);
else
    % Use the V_in from the inputs
    % disp(['getting signal from inputs'])
end
dSVdt = A * SV + B * V_in;

end