function [A] = getSS_A3(SIM,N,P,FLAG)
% Rephrase Parameters
    C_1 = SIM.C_1;
    C_4 = SIM.C_4;
    R_1 = SIM.R_1;
    R_2 = SIM.R_2;
    R_3 = SIM.R_3;
    R_4 = SIM.R_4;

    A = [-(1/R_1+1/R_2),   1/R_2       ,    0
         - 1/R_2       ,  (1/R_2+1/R_3),  - 1/R_3
          0            ,   1/R_3       ,  -(1/R_3+1/R_4)];

end