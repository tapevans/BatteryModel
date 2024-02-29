%% D_o_Liion

function [D_o_Liion] = D_oLiion_NREL(Ce , T)
    D_00 = -0.5688226;
    D_01 = -1607.003;
    D_10 = -0.8108721;
    D_11 =  475.291;
    D_20 = -0.005192312;
    D_21 = -33.43827;
    T_g0 = -24.83763;
    T_g1 =  64.07366;
    
    % D_o_Liion = (10.0.^(D_00 + D_01 ./ (T - (T_g0 + T_g1 * Ce)) + (D_10 + D_11 ./ (T - (T_g0 + T_g1 * Ce))).*Ce + (D_20 + D_21 ./ (T - (T_g0 + T_g1 * Ce))).*(Ce.^2.0)))*0.0001;
    % 
    % D_o_Liion = D_o_Liion * 4.0;

    D_o_Liion =  (D_00 + D_01 ./ (T - (T_g0 + T_g1 * Ce)))           ...
               + (D_10 + D_11 ./ (T - (T_g0 + T_g1 * Ce))).*Ce       ...
               + (D_20 + D_21 ./ (T - (T_g0 + T_g1 * Ce))).*(Ce.^2.0);
    D_o_Liion =  10.^D_o_Liion;
    D_o_Liion =  D_o_Liion * 1e-4; % Converts from cm^2/s to m^2/s
end