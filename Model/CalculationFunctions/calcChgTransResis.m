%% R_SEI Calculations
%
function [R_SEI_vec , R_SEI_vec_inv] = calcChgTransResis(CTRG , R_SEI_0)
    R_SEI_vec     = R_SEI_0 + CTRG;
    R_SEI_vec_inv = R_SEI_vec.^(-1);
end