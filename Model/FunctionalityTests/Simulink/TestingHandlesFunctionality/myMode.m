%% Testing Selecting MO Using FH
function [MO,I_ref,V_ref] = myMode(t,SIM)
fh = str2func(SIM.ProfileFnH);
[MO,I_ref,V_ref] = fh(t);
end