%% Function to calculate current flux 
% The ith index is the flux at the left (minus) face of the ith control volume
function [i_ed , i_el ] = currentCalc(  sigma_vec_interface, kappa_vec_interface, T_interface, activity_vec_interface, tf_vec_interface, ...
    phi_ed_grad, phi_el_grad, Ce_log_grad, ...
    CONS__F_inv, CONS__R, ...
    N__CV_Region_AN, N__CV_Region_SEP, N__CV_Region_CA, N__N_CV_tot, i_user)


    i_ed = -  sigma_vec_interface.*phi_ed_grad;
    i_el = -  kappa_vec_interface.*phi_el_grad ...
                    -2*kappa_vec_interface.*(CONS__R*T_interface*CONS__F_inv).*(1 + activity_vec_interface ).*(tf_vec_interface-1).*Ce_log_grad;

    i_ed(N__CV_Region_AN (1)) = i_user;
    i_el(N__CV_Region_AN (1)) = 0;
    i_ed(N__CV_Region_SEP(1)) = 0;
    i_ed(N__CV_Region_CA (1)) = 0;
    i_ed(N__N_CV_tot + 1    ) = i_user;
    i_el(N__N_CV_tot + 1    ) = 0;
end