%% Function to calculate J_Liion
% The ith index is the flux at the left (minus) face of the ith control volume
function J_Liion = JLiionCalc( D_o_Li_ion_vec_interface , tf_vec_interface, Ce_interface , ...
                               Ce_grad , T_grad , i_el , T_interface_inv_vec             , ...
                               N__CV_Region_AN, N__N_CV_tot, CONS__F_inv , EL__ST              )


    J_Liion = - D_o_Li_ion_vec_interface .* Ce_grad ...
              + tf_vec_interface * CONS__F_inv .* i_el ...
              - D_o_Li_ion_vec_interface .* Ce_interface .* EL__ST .* T_interface_inv_vec .* T_grad;

    J_Liion(N__CV_Region_AN(1)) = 0;
    J_Liion(N__N_CV_tot + 1)    = 0;
end