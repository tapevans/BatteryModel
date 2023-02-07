function [B] = getSS_B5(SIM,N,P,FLAG)
    B = zeros(5,1);
    B(P.i_PS,1) = -1;
end
