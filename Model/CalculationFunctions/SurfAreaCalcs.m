%% Calculate Surface Area Properties
%
function [A_surf_CV_vec, A_surf_CV_vec_inv , A_s_vec ] = SurfAreaCalcs(N_Particles , SIM__A_surf_single_vec , SIM__dVol_inv)
    A_surf_CV_vec     = SIM__A_surf_single_vec .* N_Particles;
    A_surf_CV_vec_inv = A_surf_CV_vec.^(-1);
    A_s_vec           = A_surf_CV_vec .* SIM__dVol_inv;
end