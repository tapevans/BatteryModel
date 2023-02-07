function [A] = getSS_A5(SIM,N,P,FLAG)
    % Rephrase Parameters
    C_1 = SIM.C_1;
    C_4 = SIM.C_4;
    R_1 = SIM.R_1;
    R_2 = SIM.R_2;
    R_3 = SIM.R_3;
    R_4 = SIM.R_4;
    
    A = zeros(5);

    A(P.V_1,P.V_1)  = - (1/R_1 + 1/R_2);
    A(P.V_1,P.V_2)  =            1/R_2 ;
    A(P.V_1,P.V_s)  =    1/R_1         ;
    
    A(P.V_2,P.V_1)  = -  1/R_2         ;
    A(P.V_2,P.V_2)  =   (1/R_2 + 1/R_3);
    A(P.V_2,P.V_3)  = -          1/R_3 ;
    
    A(P.V_3,P.V_2)  =    1/R_3         ;
    A(P.V_3,P.V_3)  = - (1/R_3 + 1/R_4);
    
    A(P.V_s,P.V_1)  =    1/R_1 ;
    A(P.V_s,P.V_s)  = -  1/R_1 ;
    A(P.V_s,P.i_PS) = -  1 ;
    
    A(P.i_PS,P.V_s) =   1 ;

end