%% Get All State Space Matricies
function [A,B,C,D] = getAll_SS5(SIM,N,P,FLAG)
    [A]      = getSS_A5(SIM,N,P,FLAG);
    [B]      = getSS_B5(SIM,N,P,FLAG);
    [~, C_m] = getSS_C5(SIM,N,P,FLAG);
               C = C_m;
    [D]      = getSS_D5(SIM,N,P,FLAG);
end