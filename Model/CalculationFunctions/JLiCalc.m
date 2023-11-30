%% Function to calculate J_Li
% The ith index is for the ith CV
% The jth index is for the flux at the left (minus) face of the jth radial control volume
% j = 1 is r = 0

function J_Li = JLiCalc( C_Li_diff, D_o_vec_interface, N__N_CV_tot, N__N_R_AN, N__N_R_CA, N__CV_Region_AN, N__CV_Region_CA, SIM__del_CV_r_inv_mat, s_dot)
    J_Li = -D_o_vec_interface .* C_Li_diff .* SIM__del_CV_r_inv_mat;

    % Symmetry BC
        j = 1;
        J_Li(j,:) = zeros( 1 , N__N_CV_tot );

    % Surface Reaction BC
        j = N__N_R_AN + 1;
        J_Li(j,N__CV_Region_AN) = s_dot(1,N__CV_Region_AN);
    
        j = N__N_R_CA + 1;
        J_Li(j,N__CV_Region_CA) = s_dot(1,N__CV_Region_CA);
end