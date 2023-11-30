%% thermalAllVectors
%
%

function [q_cond , q_conv , q_gen] = thermalAllVectors( T, T_diff, lambda_vec, N__CV_Region_AN, N__N_CV_tot, SIM__del_x_vec_halved_inv, SIM__h, SIM__T_inf, SIM__q_AN_BC, SIM__q_CA_BC, SIM__h_AN_BC, SIM__h_CA_BC, FLAG__T_BC_AN, FLAG__T_BC_CA, FLAG__HEAT_GEN_TOTAL, i_el, i_ed, i_Far)
%% Conduction 
    [q_cond] = thermalFluxConduction( T, T_diff, lambda_vec , SIM__del_x_vec_halved_inv, N__CV_Region_AN, N__N_CV_tot, SIM__T_inf, SIM__q_AN_BC, SIM__q_CA_BC, SIM__h_AN_BC, SIM__h_CA_BC, FLAG__T_BC_AN, FLAG__T_BC_CA); % [J s^-1 m^-2]


%% Convection
    [q_conv] = thermalFluxConvection( T, SIM__h, SIM__T_inf); % [J s^-1 m_outside^-2]


%% Heat Gen
    if FLAG__HEAT_GEN_TOTAL 
        % % q_gen = .0025*ones(1,N.N_CV_tot);
        % [q_gen] = thermalGenerationTotal( SV , AN , SEP, CA , EL , SIM , P , N , CONS , FLAG , props, i_el, i_ed, i_Far); % [J s^-1 ]
    else
        q_gen = zeros( 1 , N__N_CV_tot );
    end

end