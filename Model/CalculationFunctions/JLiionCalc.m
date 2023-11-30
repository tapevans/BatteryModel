%% Function to calculate J_Liion
% The ith index is the flux at the left (minus) face of the ith control volume
function J_Liion = JLiionCalc( Ce_diff, D_o_Li_ion_vec, tf_vec_interface, N__CV_Region_AN, N__N_CV_tot, CONS__F_inv, SIM__MassFluxPreCalcResistor, i_el)
    Resistor_vector = (SIM__MassFluxPreCalcResistor .* D_o_Li_ion_vec).^-1;
    Resistor_vector_tot = [nan Resistor_vector(2:end) + Resistor_vector(1:end-1)  nan];
    % J_Liion = - Resistor_vector_tot.^-1  .* Ce_diff...
    %                          + i_el .* tf_vec_interface / CONS__F;
    J_Liion = - Resistor_vector_tot.^-1  .* Ce_diff...
                             + i_el .* tf_vec_interface * CONS__F_inv;
    J_Liion(N__CV_Region_AN(1)) = 0;
    J_Liion(N__N_CV_tot + 1)    = 0;
end        