%% Get All State Space Matricies
function [A,B,C,D] = getAll_SS3(SIM,N,P,FLAG)
    [A]      = getSS_A3(SIM,N,P,FLAG);
    [B]      = getSS_B3(SIM,N,P,FLAG);
    [~, C_m] = getSS_C3(SIM,N,P,FLAG);
               C = C_m;
    [D]      = getSS_D3(SIM,N,P,FLAG);
end